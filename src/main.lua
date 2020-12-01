require "trajGen"

-------------------------------------------------------------------------------------------------------------
-- cycle time and linear and rotational trapezoidal velocity profiles
cycleTime = .001  -- s

linVel = 25    -- mm/s
linAcc = 250   -- mm/s^2
linDec = 250   -- mm/s^2

rotVel = 0.2   -- rad/s
rotAcc = 2     -- rad/s^2
rotDec = 2     -- rad/s^2

-------------------------------------------------------------------------------------------------------------
-- Define start joint position in radians
startJointPos = {0, 0, 0, 0, 0, 0}

-------------------------------------------------------------------------------------------------------------
-- Define target position, maybe joint or cartesian space
--targetPos = {0, 10, -10, rad(170), rad(30), rad(150)}  targetIsJoint = 1  -- targetPos is in joint space
targetPos = {300, 10, 300, rad(10), rad(20), rad(30)}  targetIsJoint = nil  -- targetPos is in cartesian space (x, y, z, Rx, Ry, Rz)

-------------------------------------------------------------------------------------------------------------
-- Generate trajectory, which will output joint position commands to the console
linTrajGen(startJointPos, targetPos, targetIsJoint, cycleTime, linVel, linAcc, linDec, rotVel, rotAcc, rotDec)
