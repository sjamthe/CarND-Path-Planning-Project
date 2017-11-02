# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program

### Results

Car was able to travel 5.9 miles (> 1 loop) without any incident as evident by this image and [video](https://youtu.be/QM4aEGEshes).

[![Project Video](data/final.png)](https://youtu.be/QM4aEGEshes)

### Approach
#### Summary of the solution 
function max_speed_inlane finds out how fast a car can go in the current lane (using sensor fusion data)
function change_lane decides if lane change is necessary, possible and returns new lane.
function projectPath uses the above two values and car state to determine next path for the car (using Frenet)
function smoothCoordinates uses spline to smooth the path projected above. (There is a bug in the function that I am fixing, when the car is going vertical and about to turn left (in reverse direction) spline fails as it expects sorted x values).

#### Function max_speed_inlane - How fast the car can go in current lane.
1. Use sensor fusion data to look at all the cars that are ahead in our lane.
2. For every car that is ahead of us we want to make sure we have 10 seconds gap between us and them. Based on this we determine what is the minimum speed we can take. We don't just look at car in front of us as that car may be speeding and about to hit the car in its path that has stopped!.
3. Also make our speed is not greater than the car ahead of us.
4. Adjust this speed such that each speed increment (in .02 seconds) is not very high to avoid abrupt acceleration. This limit is set as a constant max_speed_change.
Based the above we determine the speed of the car.

#### Function change_lane
This function is called for every loop is car is moving at more than 10 m/s (we don't want to change lane when car is moving slowly'). 
We use sensor fusion data to find out car ahead and behind our position in each lane and their speed.
If the our car is not moving at max speed and sees a car ahead of it that is 20 seconds away (based on current speed) we use the below logic to decide how to change lane.
1. We prefer overtaking from left, do we have a lane in our left?
2. Is the car ahead in left lane going at or above our speed?
3. Is the car in left lane atleast 20 seconds away based on our current speed (not relative speed as it may abruptly breakdown)?
4. Is the car in left lane moving faster than car in our lane?
5. Is car in left lane further ahead further away from car in front of us (may be redundant).
6. Do we have enough room in left lane to squeeze in? That means is there a car behind us in left lane? That car should be atleast 100 meters behind and should not be speeding such that it makes up that distance in 100 seconds (i.e is not faster than us by 1 m/s or more).
If all this is ok we decide to change lane to left lane.

If left lane change is not possible repeat same logic for right lane.

#### Function ProjectPath
This function takes new lane change (car_d) and new speed as input along with car coordinates and previous path.
We reuse upto 20 points from previous path, but if we are supposed to reduce speed we only use 10 points for quick reaction.
Use Frenet car coordinates and find upto 100 points ahead of the car and convert these coodinates to global coordinates.
This way we can test full move of the car from these values. The car will produce a lot of jerks that need to be smoothened.

#### Function smoothCoordinates
This function takes output of ProjectPath but only takes every 10th point of the global coordinate and creates a spline curve.
We then use this spline to find smooth global coordinates.


### Simulator.
You can download the Term3 Simulator which contains the Path Planning Project from the [releases tab (https://github.com/udacity/self-driving-car-sim/releases).

### Goals
In this project your goal is to safely navigate around a virtual highway with other traffic that is driving +-10 MPH of the 50 MPH speed limit. You will be provided the car's localization and sensor fusion data, there is also a sparse map list of waypoints around the highway. The car should try to go as close as possible to the 50 MPH speed limit, which means passing slower traffic when possible, note that other cars will try to change lanes too. The car should avoid hitting other cars at all cost as well as driving inside of the marked road lanes at all times, unless going from one lane to another. The car should be able to make one complete loop around the 6946m highway. Since the car is trying to go 50 MPH, it should take a little over 5 minutes to complete 1 loop. Also the car should not experience total acceleration over 10 m/s^2 and jerk that is greater than 50 m/s^3.

#### The map of the highway is in data/highway_map.txt
Each waypoint in the list contains  [x,y,s,dx,dy] values. x and y are the waypoint's map coordinate position, the s value is the distance along the road to get to that waypoint in meters, the dx and dy values define the unit normal vector pointing outward of the highway loop.

The highway's waypoints loop around so the frenet s value, distance along the road, goes from 0 to 6945.554.

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./path_planning`.

Here is the data provided from the Simulator to the C++ Program

#### Main car's localization Data (No Noise)

["x"] The car's x position in map coordinates

["y"] The car's y position in map coordinates

["s"] The car's s position in frenet coordinates

["d"] The car's d position in frenet coordinates

["yaw"] The car's yaw angle in the map

["speed"] The car's speed in MPH

#### Previous path data given to the Planner

//Note: Return the previous list but with processed points removed, can be a nice tool to show how far along
the path has processed since last time. 

["previous_path_x"] The previous list of x points previously given to the simulator

["previous_path_y"] The previous list of y points previously given to the simulator

#### Previous path's end s and d values 

["end_path_s"] The previous list's last point's frenet s value

["end_path_d"] The previous list's last point's frenet d value

#### Sensor Fusion Data, a list of all other car's attributes on the same side of the road. (No Noise)

["sensor_fusion"] A 2d vector of cars and then that car's [car's unique ID, car's x position in map coordinates, car's y position in map coordinates, car's x velocity in m/s, car's y velocity in m/s, car's s position in frenet coordinates, car's d position in frenet coordinates. 

## Details

1. The car uses a perfect controller and will visit every (x,y) point it recieves in the list every .02 seconds. The units for the (x,y) points are in meters and the spacing of the points determines the speed of the car. The vector going from a point to the next point in the list dictates the angle of the car. Acceleration both in the tangential and normal directions is measured along with the jerk, the rate of change of total Acceleration. The (x,y) point paths that the planner recieves should not have a total acceleration that goes over 10 m/s^2, also the jerk should not go over 50 m/s^3. (NOTE: As this is BETA, these requirements might change. Also currently jerk is over a .02 second interval, it would probably be better to average total acceleration over 1 second and measure jerk from that.

2. There will be some latency between the simulator running and the path planner returning a path, with optimized code usually its not very long maybe just 1-3 time steps. During this delay the simulator will continue using points that it was last given, because of this its a good idea to store the last points you have used so you can have a smooth transition. previous_path_x, and previous_path_y can be helpful for this transition since they show the last points given to the simulator controller with the processed points already removed. You would either return a path that extends this previous path or make sure to create a new path that has a smooth transition with this last path.

## Tips

A really helpful resource for doing this project and creating smooth trajectories was using http://kluge.in-chemnitz.de/opensource/spline/, the spline function is in a single hearder file is really easy to use.

---

## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets 
    cd uWebSockets
    git checkout e94b6e1
    ```

## Editor Settings

We've purposefully kept editor configuration files out of this repo in order to
keep it as simple and environment agnostic as possible. However, we recommend
using the following settings:

* indent using spaces
* set tab width to 2 spaces (keeps the matrices in source code aligned)

## Code Style

Please (do your best to) stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html).

## Project Instructions and Rubric

Note: regardless of the changes you make, your project must be buildable using
cmake and make!


## Call for IDE Profiles Pull Requests

Help your fellow students!

We decided to create Makefiles with cmake to keep this project as platform
agnostic as possible. Similarly, we omitted IDE profiles in order to ensure
that students don't feel pressured to use one IDE or another.

However! I'd love to help people get up and running with their IDEs of choice.
If you've created a profile for an IDE that you think other students would
appreciate, we'd love to have you add the requisite profile files and
instructions to ide_profiles/. For example if you wanted to add a VS Code
profile, you'd add:

* /ide_profiles/vscode/.vscode
* /ide_profiles/vscode/README.md

The README should explain what the profile does, how to take advantage of it,
and how to install it.

Frankly, I've never been involved in a project with multiple IDE profiles
before. I believe the best way to handle this would be to keep them out of the
repo root to avoid clutter. My expectation is that most profiles will include
instructions to copy files to a new location to get picked up by the IDE, but
that's just a guess.

One last note here: regardless of the IDE used, every submitted project must
still be compilable with cmake and make./

## How to write a README
A well written README file can enhance your project and portfolio.  Develop your abilities to create professional README files by completing [this free course](https://www.udacity.com/course/writing-readmes--ud777).

