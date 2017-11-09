#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include <algorithm>

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }
// The max s value before wrapping around the track back to 0
double MAX_S = 6945.554;
double MAX_SPEED = 49.5*1600/3600; //is 20 meters/sec
double MAX_ACC = 0.12;
double NORM_DEC = 0.12;
double MAX_DEC = 0.22;
double SAFE_GAP = 20;

//double max_speed_change = 0.12; //don't change speed too much to avoid max out accelerations

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0; //SJ: was 0

    double dist;
	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
        dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}
	}

	return closestWaypoint;

}


int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
    
    int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);
    int nextWaypoint;
    
    double wp_x = maps_x[closestWaypoint];
    double wp_y = maps_y[closestWaypoint];
    
    double prev_x, prev_y;
    if(closestWaypoint > 0) {
         prev_x = maps_x[closestWaypoint-1];
         prev_y = maps_y[closestWaypoint-1];
    } else {
        //circular track
        prev_x = maps_x[maps_x.size()-1];
        prev_y = maps_y[maps_x.size()-1];
    }
    
    double prev_wpdist, next_wpdist, dist1, dist2, dist3;
    
    prev_wpdist = distance(wp_x, wp_y, prev_x, prev_y);
    //next_wpdist = distance(wp_x, wp_y, next_x, next_y);
    dist1 = distance(x, y, prev_x, prev_y);
    dist2 = distance(wp_x, wp_y, x, y);

    //If the angle is obtuse then car is ahead of the waypoint.
    if(pow(dist1,2) - (pow(dist2,2)  + pow(prev_wpdist,2)) > 0.1 || //We use 0.1 to avoid rounding 0 is the right value
       pow(dist2,2) - (pow(dist1,2)  + pow(prev_wpdist,2)) > 0.1) {
        nextWaypoint = closestWaypoint+1;
    } else {
        nextWaypoint = closestWaypoint;
    }
    if(nextWaypoint == maps_x.size()) {
        nextWaypoint = 0;
    }
    //cout << "returning next wp " << nextWaypoint << endl;
    
    return nextWaypoint;
    
}
// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x, maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}
    //cout << "getFrenet: prev_wp = " << prev_wp << " , next_wp = " <<  next_wp  << " for x, y, theta " << x << ", " << y << ", " << theta << endl;
	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value -- SJ : Shouldn't we use map_s instead of calculating distance from 0?
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

    //Find the way point that is just in front of the car (s)
	while(s > maps_s[prev_wp+1] && (prev_wp + 1 < (int)(maps_s.size()) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
    //cout << "getXY: waypoints used " << prev_wp <<"," << wp2 << " heading = " << heading ;
    //cout << " car at " << s << " between " << maps_s[prev_wp] << "," << maps_s[wp2] << endl;
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

    double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

//Convert global coordinates to car-coordinates
vector<double> carCoordinates(double x, double y, double ref_x, double ref_y, double ref_yaw) {
    //translation of points
    double dx = x - ref_x;
    double dy = y - ref_y;
    //rotation
    x = dx * cos(ref_yaw) + dy * sin(ref_yaw);
    y = dy * cos(ref_yaw) - dx * sin(ref_yaw);
    /*if(x < 0) {
        cout << " firection " << x << endl;
    }*/
    return {x, y};
}

//Convert car coordinates to global coordinates
vector<double> globalCoordinates(double x_point, double y_point, double ref_x, double ref_y, double ref_yaw) {
    //convert vehicle coordinates back to regular system x,y coordinates
    double x = x_point*cos(ref_yaw) - y_point*sin(ref_yaw);
    double y = x_point*sin(ref_yaw) + y_point*cos(ref_yaw);
    //now translate it back.
    x += ref_x;
    y += ref_y;
    
    return {x, y};
}

//sort x,y vectors by x as spline needs sorted x
vector<vector<double>>  sortXY(vector <double> x_vals, vector < double> y_vals) {
    vector <pair <double, double> > pairxy;
    for (int i=0; i<x_vals.size(); i++) {
        pairxy.push_back( make_pair(x_vals[i], y_vals[i]) );
    }
    sort(pairxy.begin(), pairxy.end());
    
    vector<double> x;
    vector<double> y;
    for (int i=0; i<x_vals.size(); i++) {
      if(x.size() == 0 || pairxy[i].first > x[x.size()-1]) {
        x.push_back(pairxy[i].first);
        y.push_back(pairxy[i].second);
      }
    }
    return {x, y};
}

// projectPath: gets new lane change (car_d is needed) and new speed as input along with car coordinates and previous values.
// After taking some previous values for smooth transition it predict new path.
// We first deternime the cars next path in Frenet. We increase car_s by velocity*0.02 and also adjust car_d slightly over 100 points
// The frenet coordinates are conveted to global coordinates and returned.
vector<vector<double>> projectPath(double car_x, double car_y, double car_s, double car_d, double new_d, double car_yaw,
                                   double car_speed, double final_speed,
             vector<double> previous_path_x, vector<double> previous_path_y, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
    vector<double> x_vals;
    vector<double> y_vals;
    
    int total_points = 50;
    
    double prev_size = previous_path_x.size();
    
    //cout << "prev_size = " << prev_size << endl;
    //initialize ref to car's location. this changes only if we use prev points
    double ref_x = car_x;
    double ref_y = car_y;
    double ref_yaw = car_yaw;
    
    //start new path aligned with the angle of the car
    if(prev_size < 2) {
        double prev_car_x = car_x - cos(car_yaw);
        double prev_car_y = car_y - sin(car_yaw);
        
        x_vals.push_back(prev_car_x);
        x_vals.push_back(car_x);
        
        y_vals.push_back(prev_car_y);
        y_vals.push_back(car_y);
    }
    else { //start with two last points from prev_path
        ref_x = previous_path_x[prev_size-1];
        ref_y = previous_path_y[prev_size-1];
        
        double prev_ref_x = previous_path_x[prev_size-2];
        double prev_ref_y = previous_path_y[prev_size-2];
        ref_yaw = atan2(ref_y-prev_ref_y, ref_x-prev_ref_x);
        
        x_vals.push_back(prev_ref_x);
        x_vals.push_back(ref_x);
        
        y_vals.push_back(prev_ref_y);
        y_vals.push_back(ref_y);
    }
    vector<double> sd = getFrenet(ref_x, ref_y, ref_yaw, maps_s, maps_x, maps_y);
    cout << "sd[0], car_s " << sd[0] << ", " << car_s << endl;
    
    //Add 30m spaced points ahead of the car (shouldn't this be ref_s?)
    vector<double> next_xy0 = getXY(sd[0]+30, new_d, maps_s, maps_x, maps_y);
    vector<double> next_xy1 = getXY(sd[0]+60, new_d, maps_s, maps_x, maps_y);
    vector<double> next_xy2 = getXY(sd[0]+90, new_d, maps_s, maps_x, maps_y);
    
    x_vals.push_back(next_xy0[0]);
    x_vals.push_back(next_xy1[0]);
    x_vals.push_back(next_xy2[0]);
    
    y_vals.push_back(next_xy0[1]);
    y_vals.push_back(next_xy1[1]);
    y_vals.push_back(next_xy2[1]);
    
    //Transform coordinates to local car coordinates
    for(int i=0; i<x_vals.size(); i++) {
        //shift car ref angle to zero degree
        double shift_x = x_vals[i] - ref_x;
        double shift_y = y_vals[i] - ref_y;
        
        //rotate the point
        x_vals[i] = (shift_x*cos(0-ref_yaw) - shift_y*sin(0-ref_yaw));
        y_vals[i] = (shift_x*sin(0-ref_yaw) + shift_y*cos(0-ref_yaw));
    }
    
    tk::spline s;
    s.set_points(x_vals,y_vals,true); //true is for cubic spline
    
    //predict car waypoints
    vector<double> next_x_vals;
    vector<double> next_y_vals;
    
    //Add max 20 previous_path points
    //int start = max(prev_size-20.0, 0.0);
    for (int i=0; i<prev_size; i++) {
        next_x_vals.push_back(previous_path_x[i]);
        next_y_vals.push_back(previous_path_y[i]);
    }
    
    //get points from spline
    //let us assume that 30m on spline is almost straight.
    double target_x = 30;
    double target_y = s(target_x);
    double target_s = sqrt(target_x*target_x + target_y*target_y);
    double N = target_s/(.02*final_speed); //Number of waypoints in this distance.
    
    //we start taking points from begin of spline which started at car_s (but if we have prev points that should fail)
    //cout << "N = " << N  << endl;
    double current_x = 0;
    for(int i=1; i<=total_points-prev_size; i++) {
        double incr_x = 0.02*final_speed;
        double local_x = current_x + incr_x;
        double local_y = s(local_x);
        
        //increment values for next loop
        current_x += incr_x;

        //rotate the point back
        double x_point = (local_x*cos(ref_yaw) - local_y*sin(ref_yaw));
        double y_point = (local_x*sin(ref_yaw) + local_y*cos(ref_yaw));
        
        //shift point back
        x_point += ref_x;
        y_point += ref_y;
        
        next_x_vals.push_back(x_point);
        next_y_vals.push_back(y_point);
    }
    
    return {next_x_vals, next_y_vals};
}

//Use sensor_fusion data to find if it is safe to change lane.
//we are looking for a lane gap. relative to our car_s find out car ahead and behind in each lane.
//
double change_lane(vector<vector <double>> sensor_fusion, double car_s, double car_d, double car_speed, double final_speed) {
    
    int car_lane;
    
    //use a little save zones for lanes
    if(car_d <= 4.5) {
        car_lane = 0;
    } else if(car_d > 3.5 && car_d <= 8.5) {
        car_lane = 1;
    } else {
        car_lane = 2;
    }
    
    //stores location of car that is closest (but in front) to car_s in each lane.
    vector <double> closet_car_ahead = {7000, 7000, 7000};
    //speed of the car in front
    vector <double> speed_car_ahead = {MAX_SPEED, MAX_SPEED, MAX_SPEED};
    
    //stores location of car that is closest (but in back) to car_s in each lane.
    vector <double> closet_car_behind = {7000, 7000, 7000};
    //speed of the car behind
    vector <double> speed_car_behind = {MAX_SPEED, MAX_SPEED, MAX_SPEED};
    
    for(int i = 0; i < sensor_fusion.size();i++) {
        // car is in my lane
        int sf_id = sensor_fusion[i][0];
        double sf_x = sensor_fusion[i][1];
        double sf_y = sensor_fusion[i][2];
        double sf_vx = sensor_fusion[i][3];
        double sf_vy = sensor_fusion[i][4];
        double sf_s = sensor_fusion[i][5];
        double sf_d = sensor_fusion[i][6];
        double sf_speed = sqrt(sf_vx * sf_vx + sf_vy * sf_vy); //speed of the car.
        double sf_car_dist = sf_s - car_s; //distance of this car from our car.
        
        //handle circular track situation
        if(sf_car_dist > MAX_S/2) {
            //we are crossing the loop but sf_car is actually behind us.
            cout << sf_d << " before adjustment " << sf_car_dist << " after " ;
            sf_car_dist = MAX_S - sf_s + car_s;
            cout << sf_car_dist << endl;
        } else if (sf_car_dist < -1*MAX_S/2) {
            //we are behind the 0 point but sf_s is ahead of us.
            cout << sf_d << " before adjustment " << sf_car_dist << " after " ;
            sf_car_dist = MAX_S - car_s + sf_s;
            cout << sf_car_dist << endl;
        }
        
        int sf_lane = int(sf_d/4); //lane of this car
        
        if(sf_car_dist > 0) { //car is ahead
            //A car can be in more than one lane if it is changing lane. be careful of these cars
            //and add them in both lanes.
            for(int lane=0;lane<3;lane++) {
                double d_min = lane*4+2-3;//instead of +=2 we do +=3 to cover lane changers
                double d_max = lane*4+2+3;
                //cout << "ahead in lane " << lane << " at " << sf_car_dist << d_min << ", " << d_max << endl;
                if(sf_d >= d_min && sf_d <= d_max) {
                    if(sf_car_dist < closet_car_ahead[lane]) {
                        closet_car_ahead[lane] = sf_car_dist;
                        speed_car_ahead[lane] = sf_speed;
                    }
                }
            }
        } else { //car is behind
            cout << "car behind at " << abs(sf_car_dist) << " , lane = " << sf_lane << endl;
            if(abs(sf_car_dist) < closet_car_behind[sf_lane]) {
                closet_car_behind[sf_lane] = abs(sf_car_dist);
                speed_car_behind[sf_lane] = sf_speed;
            }
        }
    }
    //debug
    for(int i=0; i<3; i++) {
        cout << "lane " << i << " car ahead dist speed, " << closet_car_ahead[i] << ", " << speed_car_ahead[i] <<
        ", car behind dist speed, "<< closet_car_behind[i] << ", " << speed_car_behind[i] << endl;
    }
    
    //Is there a need to change lane? if speed is not max and we are blocked ahead
    //If we need to change lane, prefer left lane over right.
    //make sure car in that lane is moving faster than us.
    //and there is enough gap (30 seconds?) in front and back to change lane.
    double prev_lane = car_lane;
    int front_gap = 35;
    int rear_gap = 30;
    int prefer_lane = 1;
    
    if((final_speed <= MAX_SPEED && closet_car_ahead[car_lane] < front_gap) || car_lane != prefer_lane) {
        cout << "slow speed and car in front " << closet_car_ahead[car_lane] << " > front_gap , ";
        if(car_lane > 0) { //we can try left lane first
            cout << " left, " ;
            if(speed_car_ahead[car_lane-1] > speed_car_ahead[car_lane]) {
                cout << " fast, ";
                //see if front gap in that lane is greater that current lane
                if(closet_car_ahead[car_lane-1] > front_gap ) {
                    cout << " room ahead, ";
                    if(closet_car_behind[car_lane - 1] > rear_gap) {
                        if ((closet_car_behind[car_lane - 1] >= MAX_S/2) || //car too far
                            (speed_car_behind[car_lane - 1] - final_speed)*5 < rear_gap) { //car won't come close
                           cout << " change lane to left " << closet_car_behind[car_lane - 1] << " > " << rear_gap <<
                           " && speed diff " << (speed_car_behind[car_lane - 1] - final_speed);
                           car_lane = car_lane - 1;
                       } else {
                           cout << car_lane << " don't change to left as car behind too fast " << (speed_car_behind[car_lane - 1] - final_speed)*10 <<
                           " > " << rear_gap ;
                       }
                    } else {
                        cout << "don't change to left as car behind too close " << closet_car_behind[car_lane - 1] << " < " << rear_gap;
                    }
                } else {
                    cout << "car ahead on left lane close at " << closet_car_ahead[car_lane-1] << " < " << front_gap ;
                }
            } else {
                cout << " car in left lane at speed " << speed_car_ahead[car_lane-1] << " < " << speed_car_ahead[car_lane] ;
            }
        }
        //if we couldn't change to left, try to right.
        if(prev_lane == car_lane && car_lane < 2) {
            cout << " right, " ;
            if(speed_car_ahead[car_lane+1] > speed_car_ahead[car_lane]) {
                cout << " fast, ";
                //see if front gap in that lane is greater that current lane
                if(closet_car_ahead[car_lane+1] > front_gap) {
                    cout << " room ahead, ";
                    if(closet_car_behind[car_lane + 1] > rear_gap) {
                        if((closet_car_behind[car_lane + 1] >= MAX_S/2) || //car too far
                           (speed_car_behind[car_lane + 1] - final_speed)*5 < rear_gap) { //car won't come close
                            cout << " change lane to right " << closet_car_behind[car_lane + 1] << " > " << rear_gap <<
                            " && speed diff " << (speed_car_behind[car_lane + 1] - final_speed);
                            car_lane = car_lane + 1;
                        } else {
                            cout << "don't change to right as car behind too fast " << (speed_car_behind[car_lane + 1] - final_speed)*10 <<
                            " > " << rear_gap ;
                        }
                    } else {
                        cout << "don't change to right as car behind too close " << closet_car_behind[car_lane + 1] << " < " << rear_gap;
                    }
                } else {
                    cout <<  "car ahead on right lane at " << closet_car_ahead[car_lane+1] << " < " << front_gap << endl;
                }
            } else {
                cout << "slow car in right lane at " << speed_car_ahead[car_lane+1] << " at speed " << speed_car_ahead[car_lane];
            }
        }
    } else {
        cout << "no need to change lane as car ahead at " << closet_car_ahead[car_lane] << " > " <<  front_gap << endl;
    }
    car_d = car_lane*4 + 2;
    cout << " prev_lane, car_lane " << prev_lane << ", " << car_lane << " car_d " << car_d << endl;
    return car_d;
}

//Use sensor_fusion data to find out what is the max speed we can move in this lane.
double max_speed_inlane(vector<vector <double>> sensor_fusion, double car_s, double car_d, double car_speed) {
  
    double final_speed = car_speed;
    bool slowdown = false;
    double closest_car_dist = 7000;
    double closest_car_speed = MAX_SPEED;

    cout << "start final = " << final_speed;
    for(int i = 0; i < sensor_fusion.size();i++) {
        // car is in my lane
        int sf_id = sensor_fusion[i][0];
        double sf_x = sensor_fusion[i][1];
        double sf_y = sensor_fusion[i][2];
        double sf_vx = sensor_fusion[i][3];
        double sf_vy = sensor_fusion[i][4];
        double sf_s = sensor_fusion[i][5];
        double sf_d = sensor_fusion[i][6];

        // see if car at same lane (lane width is 4
        if (sf_d >= car_d - 2 && sf_d <= car_d + 2 ) {
            
            if(sf_s < car_s && (car_s - sf_s) < MAX_S/2) {
                //this car is behind us.
                continue;
            }
          
            double sf_speed = sqrt(sf_vx * sf_vx + sf_vy * sf_vy);
            
            double sf_car_dist = sf_s - car_s;
            //handle circular track situation
            if(sf_car_dist < 0) {
                sf_car_dist = MAX_S - car_s + sf_s;
                cout << "handling circular track in max_speed_inlane sf_car_dist, " <<  sf_car_dist << ", "
                        << sf_s << ", " << car_s << endl;
            }
            
            if(closest_car_dist > sf_car_dist) {
                closest_car_dist = sf_car_dist;
                closest_car_speed = sf_speed;
            }
          
          //We want to keep safe distance from all cars in our lane that are ahead of us
          //that means we should have a gap such that it takes 5 seconds to cover based on relative velocity.
          //We should look at all cars not just the closet car as there may be accident ahead and closet car may not have reacted to it.
          
          //speed needed to cover the gap in 5 seconds
          double speed = (sf_car_dist - SAFE_GAP)/5 + sf_speed;
          if(final_speed > speed) {
            slowdown = true;
            final_speed = speed;
            cout << ", safegap speed " << final_speed;
          }
        }
    }
    cout << ", closest_car_speed = " << closest_car_speed << ", final speed " << final_speed ;
    //are we going faster than closest car?
    if(final_speed > closest_car_speed && closest_car_dist < SAFE_GAP) { //30 meters gap seems good
        final_speed = closest_car_speed ; //maintain speed of the car ahead.
        slowdown = true;
        cout << ", slowing down for closest car " << final_speed ;
    }
    if(slowdown) {
        //see if we have slowed down more than permitted DEC
        if(final_speed-car_speed > MAX_DEC) {
            final_speed = car_speed - MAX_DEC;
        }
    } else {
        //accelerate if we are not at the max
        if(car_speed + MAX_ACC < MAX_SPEED) {
            final_speed = car_speed + MAX_ACC;
        }
    }
    
    cout << " final = " << final_speed << endl;
    return final_speed;
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;
  double old_yaw = 0;
  double old_car_s = 0;
  double prev_change_s = 0;
  double old_speed = 0;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&old_yaw, &old_car_s,
               &old_speed, &prev_change_s](uWS::WebSocket<uWS::SERVER> ws,
                char *data, size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"]; //even yaw appears to be wrong sometimes?
          	double car_speed = j[1]["speed"]; //speed returned is not acurate, don't use it.
            car_speed = car_speed*1600/3600; //in m/s
            
            if(abs(car_speed-old_speed) > 10) {
                cout << "sudden speed change from " << old_speed << " to " << car_speed << endl;
                //car_speed = old_speed;
            }

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
            int prev_size = previous_path_x.size();
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];
            
          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;
            
            //We seem to be getting yaw very large, make it less than pi
            double ref_yaw = deg2rad(car_yaw);
          
            while (ref_yaw >= 2*pi()) {
                cout << "reducing ref_yaw by 2pi " << ref_yaw << endl;
                ref_yaw -= 2*pi();
            }
           
            /*if(car_s > MAX_S)
                car_s -= MAX_S;
             */
            //Use sensor_fusion data to find out what is the max speed we can move in this lane.
            double max_speed = max_speed_inlane(sensor_fusion, car_s, car_d, old_speed);
            
            //cout << "car_s, car_d, ref_yaw, car_speed, speed-diff, " << car_s << ", " << car_d << ", " << ref_yaw << ", " << car_speed << ", " << max_speed-car_speed << endl;
            //assert(car_s  >=old_car_s -0.75);
            
            //Allow lane change only when moving above 8 m/s
            double new_d = car_d;
            if(max_speed > 8) {
                //change lane if necessary else stick close to center of the lane car is already on.
                new_d = change_lane(sensor_fusion, car_s, car_d, car_speed, max_speed);
                prev_change_s = car_s;
            }
            /*don't change speed if changing lane to reduce accl limit
            if(int(new_d/4) != int(car_d/4)) {
                max_speed = car_speed;
                cout << "no change in speed for lane change max_speed = " << max_speed << endl;
            }*/
            
            vector<vector<double>> results = projectPath(car_x, car_y, car_s, car_d, new_d, ref_yaw, car_speed, max_speed, previous_path_x,
                                                        previous_path_y, map_waypoints_s, map_waypoints_x, map_waypoints_y);

            /* results = smoothCoordinates(results[0], results[1], car_x, car_y, ref_yaw, max_speed,
                                                        car_s, map_waypoints_s, map_waypoints_x, map_waypoints_y);*/
            msgJson["next_x"] = results[0];
            msgJson["next_y"] = results[1];
            
            old_yaw = ref_yaw;
            old_car_s = car_s;
            old_speed = max_speed;
            
          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
    
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
