require "kinematics"

-------------------------------------------------------------------------------------------------------------
--   Trapezoidal Velocity profile types
--
--    ^ V     _________________                                   ^ V                     
--    |     /|                 |\                                 |                    
--    |    / |                 | \                                |    /|\              
--    |acc/  |                 |  \dec                            |acc/ | \dec          
--    |  /   |h                |h  \                   or         |  /  |h \            
--    | /    |                 |    \                             | /   |   \           
--    |/     |                 |     \                            |/    |    \          
--    -------------------------------------> t                    ------------------> t   
--       t1          t3           t2                                 t1    t2            
--                                                                                        
--                 case 1                                             case 2
-------------------------------------------------------------------------------------------------------------
-- T  = simulation time in seconds
-- t1 = velocity rise time
-- t2 = velocity fall time
-- T = t1 + t2
--
-- Define:
-- h      = max achievable velocity during trapezoidal velocity move
-- maxVel = max achievable velocity with given position target and profile acceleration and deceleration
-- 
-- if maxVel > vel:
--   h = vel         (case 1)
-- else:
--   h = maxVel      (case 2)
--
-------------------------------------------------------------------------------------------------------------
-- Solving for maxVel
--
--    ^ V
--    |       /|\
--    |      / | \            
--    |     /  |  \            
--    |    /   |   \          
--    |acc/  maxVel \ dec     
--    |  /     |     \         
--    | /      |      \        
--    |/       |       \       
--    --------------------> t  
--        t1       t2                  
--
--                                                                                        
--         maxVel                      maxVel             
-- acc =  --------      --->    t1 =  --------
--           t1                          acc               
--
--
--         maxVel                      maxVel             
-- dec =  --------      --->    t2 =  --------
--           t2                         dec    
--
--
--         maxVel x T                 2 x dis             
-- dis =  ------------  --->    T =  ---------
--              2                      maxVel       
--                                                                                       _____________________
--                               2 x Pos     maxVel    maxVel                      \    / 2 x dis x acc x dec
-- T = t1 + t2          --->    --------- = -------  + -------     --->   maxVel =  \  /  --------------------
--                                maxVel      acc        dec                         \/       acc + dec
--

-------------------------------------------------------------------------------------------------------------
-- 3 dimensional trajectory generator

function trajGenXYZInit(pi, pf, vel, acc, dec)

  trajGenXYZ = {}

  local dx = pf[1]-pi[1]
  local dy = pf[2]-pi[2]
  local dz = pf[3]-pi[3]

  local dis = sqrt(dx*dx + dy*dy + dz*dz)

  -- 1D parameters
  local h = min(vel, sqrt(2*dis*acc*dec/(acc+dec)))
  local t1 = h/acc
  local t2 = h/dec
  local t3 = 0
  if h ~= 0 then
    t3 = (dis-(t1+t2)*h/2)/h
  end
  local t13 = t1+t3
  local T = t1+t2+t3

  trajGenXYZ["dis"] = dis
  trajGenXYZ["vel"] = vel
  trajGenXYZ["acc"] = acc
  trajGenXYZ["dec"] = dec

  trajGenXYZ["h"]   = h
  trajGenXYZ["t1"]  = t1
  trajGenXYZ["t3"]  = t3
  trajGenXYZ["t13"] = t13
  trajGenXYZ["T"]   = T


  -- 3D parameters
  trajGenXYZ["xi"] = pi[1]
  trajGenXYZ["yi"] = pi[2]
  trajGenXYZ["zi"] = pi[3]

  if dis ~= 0 then
    trajGenXYZ["kx"] = dx/dis
    trajGenXYZ["ky"] = dy/dis
    trajGenXYZ["kz"] = dz/dis
  else
    trajGenXYZ["kx"] = 0
    trajGenXYZ["ky"] = 0
    trajGenXYZ["kz"] = 0
  end
end

function trajGenXYZGetP(t, xyz)

  local p = 0
  local v = 0
  local a = 0

  if t<trajGenXYZ.t1 then
    a = trajGenXYZ.acc
    v = trajGenXYZ.acc*t
    p = 1/2*trajGenXYZ.acc*t*t
  elseif t<trajGenXYZ.t13 then
    a = 0
    v = trajGenXYZ.h
    p = 1/2*trajGenXYZ.acc*trajGenXYZ.t1*trajGenXYZ.t1+(t-trajGenXYZ.t1)*trajGenXYZ.h
  elseif t<trajGenXYZ.T then
    a = -trajGenXYZ.dec
    v = trajGenXYZ.h-trajGenXYZ.dec*(t-trajGenXYZ.t13)
    p = 1/2*trajGenXYZ.t1*trajGenXYZ.h+trajGenXYZ.t3*trajGenXYZ.h+(t-trajGenXYZ.t13)*trajGenXYZ.h-1/2*trajGenXYZ.dec*(t-trajGenXYZ.t13)*(t-trajGenXYZ.t13)
  else
    a = 0
    v = 0
    p = trajGenXYZ.dis
  end

  local px = p*trajGenXYZ.kx + trajGenXYZ.xi
  local py = p*trajGenXYZ.ky + trajGenXYZ.yi
  local pz = p*trajGenXYZ.kz + trajGenXYZ.zi

  local vx = v*trajGenXYZ.kx
  local vy = v*trajGenXYZ.ky
  local vz = v*trajGenXYZ.kz

  local ax = a*trajGenXYZ.kx
  local ay = a*trajGenXYZ.ky
  local az = a*trajGenXYZ.kz

  xyz[1] = px
  xyz[2] = py
  xyz[3] = pz
end

-------------------------------------------------------------------------------------------------------------
-- 4 dimensional trajectory generator, producing a trapezoidal velocity profile on quaternions using SLERP
-- For more details about SLERP algorithm please see: https://en.wikipedia.org/wiki/Slerp

function trajGenQUATInit(qi, qf, vel, acc, dec)

  trajGenQUAT = {}

  local dot = (qi[1]*qf[1] + qi[2]*qf[2] + qi[3]*qf[3] + qi[4]*qf[4])

  if dot < 0 then
    qf[1] = -qf[1]
    qf[2] = -qf[2]
    qf[3] = -qf[3]
    qf[4] = -qf[4]
    dot = -dot
  end

  -- Angle of rotation in quaternion space
  local omega = acos(dot)

  -- Angle of rotation in 3D space
  local dis  = 2*omega   --IMPORTANT NOTE! A rotation of omega in quaternion space maps to a rotation of 2*omega in 3D space

    -- 1D parameters
  local h = min(vel, sqrt(2*dis*acc*dec/(acc+dec)))
  local t1 = h/acc
  local t2 = h/dec
  local t3 = 0
  if h ~= 0 then
    t3 = (dis-(t1+t2)*h/2)/h
  end
  local t13 = t1+t3
  local T = t1+t2+t3

  trajGenQUAT["dis"] = dis
  trajGenQUAT["vel"] = vel
  trajGenQUAT["acc"] = acc
  trajGenQUAT["dec"] = dec

  trajGenQUAT["h"]   = h
  trajGenQUAT["t1"]  = t1
  trajGenQUAT["t3"]  = t3
  trajGenQUAT["t13"] = t13
  trajGenQUAT["T"]   = T


  -- 4D parameters
  trajGenQUAT["omega"] = omega
  trajGenQUAT["somega"] = sin(omega)

  trajGenQUAT["dot"] = dot

  trajGenQUAT["qi"] = qi
  trajGenQUAT["qf"] = qf
end

function trajGenQUATGetQ(t, q)

  local p = 0
  local v = 0
  local a = 0

  if t<trajGenQUAT.t1 then
    a = trajGenQUAT.acc
    v = trajGenQUAT.acc*t
    p = 1/2*trajGenQUAT.acc*t*t
  elseif t<trajGenQUAT.t13 then
    a = 0
    v = trajGenQUAT.h
    p = 1/2*trajGenQUAT.acc*trajGenQUAT.t1*trajGenQUAT.t1+(t-trajGenQUAT.t1)*trajGenQUAT.h
  elseif t<trajGenQUAT.T then
    a = -trajGenQUAT.dec
    v = trajGenQUAT.h-trajGenQUAT.dec*(t-trajGenQUAT.t13)
    p = 1/2*trajGenQUAT.t1*trajGenQUAT.h+trajGenQUAT.t3*trajGenQUAT.h+(t-trajGenQUAT.t13)*trajGenQUAT.h-1/2*trajGenQUAT.dec*(t-trajGenQUAT.t13)*(t-trajGenQUAT.t13)
  else
    a = 0
    v = 0
    p = trajGenQUAT.dis
  end

  SIN_THRESHOLD = 0.0005
  
  local theta = p/2 -- IMPORTANT NOTE! A rotation of p in 3D space maps to a rotation of p/2 in quaternion space
  local s1 = 0
  local s0 = 0

  if (trajGenQUAT.somega > SIN_THRESHOLD) then -- use slerp
    s1 = sin(theta) / trajGenQUAT.somega
    s0 = cos(theta) - trajGenQUAT.dot * s1
  else                                         -- use lerp
    s1 = theta/trajGenQUAT.omega
    s0 = 1 - s1
  end

  q[1] = s0*trajGenQUAT.qi[1] + s1*trajGenQUAT.qf[1]
  q[2] = s0*trajGenQUAT.qi[2] + s1*trajGenQUAT.qf[2]
  q[3] = s0*trajGenQUAT.qi[3] + s1*trajGenQUAT.qf[3]
  q[4] = s0*trajGenQUAT.qi[4] + s1*trajGenQUAT.qf[4]

  normalizeQuat(q)
end

-------------------------------------------------------------------------------------------------------------
-- 3 dimensional trajectory generator

function trajGenRBTInit(pi, pf, vel, acc, dec)

  trajGenRBT = {}

  local dx = pf[1]-pi[1]
  local dy = pf[2]-pi[2]
  local dz = pf[3]-pi[3]

  local dis = sqrt(dx*dx + dy*dy + dz*dz)

  -- 1D parameters
  local h = min(vel, sqrt(2*dis*acc*dec/(acc+dec)))
  local t1 = h/acc
  local t2 = h/dec
  local t3 = 0
  if h ~= 0 then
    t3 = (dis-(t1+t2)*h/2)/h
  end
  local t13 = t1+t3
  local T = t1+t2+t3

  trajGenRBT["dis"] = dis
  trajGenRBT["vel"] = vel
  trajGenRBT["acc"] = acc
  trajGenRBT["dec"] = dec

  trajGenRBT["h"]   = h
  trajGenRBT["t1"]  = t1
  trajGenRBT["t3"]  = t3
  trajGenRBT["t13"] = t13
  trajGenRBT["T"]   = T


  -- 3D parameters
  trajGenRBT["xi"] = pi[1]
  trajGenRBT["yi"] = pi[2]
  trajGenRBT["zi"] = pi[3]

  if dis ~= 0 then
    trajGenRBT["kx"] = dx/dis
    trajGenRBT["ky"] = dy/dis
    trajGenRBT["kz"] = dz/dis
  else
    trajGenRBT["kx"] = 0
    trajGenRBT["ky"] = 0
    trajGenRBT["kz"] = 0
  end
end

function trajGenRBTGetP(t, rbt)

  local p = 0
  local v = 0
  local a = 0

  if t<trajGenRBT.t1 then
    a = trajGenRBT.acc
    v = trajGenRBT.acc*t
    p = 1/2*trajGenRBT.acc*t*t
  elseif t<trajGenRBT.t13 then
    a = 0
    v = trajGenRBT.h
    p = 1/2*trajGenRBT.acc*trajGenRBT.t1*trajGenRBT.t1+(t-trajGenRBT.t1)*trajGenRBT.h
  elseif t<trajGenRBT.T then
    a = -trajGenRBT.dec
    v = trajGenRBT.h-trajGenRBT.dec*(t-trajGenRBT.t13)
    p = 1/2*trajGenRBT.t1*trajGenRBT.h+trajGenRBT.t3*trajGenRBT.h+(t-trajGenRBT.t13)*trajGenRBT.h-1/2*trajGenRBT.dec*(t-trajGenRBT.t13)*(t-trajGenRBT.t13)
  else
    a = 0
    v = 0
    p = trajGenRBT.dis
  end

  local px = p*trajGenRBT.kx + trajGenRBT.xi
  local py = p*trajGenRBT.ky + trajGenRBT.yi
  local pz = p*trajGenRBT.kz + trajGenRBT.zi

  local vx = v*trajGenRBT.kx
  local vy = v*trajGenRBT.ky
  local vz = v*trajGenRBT.kz

  local ax = a*trajGenRBT.kx
  local ay = a*trajGenRBT.ky
  local az = a*trajGenRBT.kz

  rbt[1] = px
  rbt[2] = py
  rbt[3] = pz
end

-------------------------------------------------------------------------------------------------------------

function linTrajGen(startJointPos, targetPos, targetIsJoint, cycleTime, linVel, linAcc, linDec, rotVel, rotAcc, rotDec)

  -- Copy start joint position
  vectorCopy(startJointPos, 1, theta1, 1, 6)

  local rslti = forward(theta1, world1)

  if rslti~=success then
    print('Error in forward kinematics at start position: ', rslti)
    exit()
  end

  if targetIsJoint then
    vectorCopy(targetPos, 1, theta2, 1, 6)

    local rsltf = forward(theta2, world2)

    if rsltf~=success then
      print('Error in forward kinematics at end position: ', rsltf)
      exit()
    end
  else
    vectorCopy(startJointPos, 1, theta2, 1, 6)

    vectorCopy(targetPos, 1, world2, 1, 3)
    vectorCopy(targetPos, 4, abg,    1, 3)
    eulerToQuat(abg, quat)

    vectorCopy(quat, 1, world2, 4, 4)

    local rsltf = inverse(world2, theta2)

    if rsltf~=success then
      print('Error in inverse kinematics at end position: ', rsltf)
      exit()
    end
  end

  if not targetIsJoint then
    chooseCloserWristConfig(theta1, theta2)      -- Flip wrist config of theta2 if it provides a shorter joint path
    constainRJointInWristSing(theta1, theta2)    -- Constrain R joint from changing if end config is in wrist singularity
  end

  local interpRBT = isWristSingInPath(theta1, theta2) -- If wrist singularity is along the path, then interpolate wrist joints

  -- xyz interpolation
  vectorCopy(world1, 1, xyz1, 1, 3)
  vectorCopy(world2, 1, xyz2, 1, 3)
  trajGenXYZInit(xyz1, xyz2, linVel, linAcc, linDec)

  local n = trajGenXYZ.T/cycleTime

  if interpRBT then -- rbt interpolation
    vectorCopy(theta1, 4, rbt1, 1, 3)
    vectorCopy(theta2, 4, rbt2, 1, 3)

    trajGenRBTInit(rbt1, rbt2, rotVel, rotAcc, rotDec)
    n = max(n, trajGenRBT.T/cycleTime)  -- Choose the longer of the two moves
  else -- spherical quat interpolation (SLERP)
    vectorCopy(world1, 4, quat1, 1, 4)
    vectorCopy(world2, 4, quat2, 1, 4)

    trajGenQUATInit(quat1, quat2, rotVel, rotAcc, rotDec)
    n = max(n, trajGenQUAT.T/cycleTime)  -- Choose the longer of the two moves
  end

  -- set initial hint position to know initial configuration
  vectorCopy(startJointPos, 1, theta, 1, 6)

  for i=1,n do
    trajGenXYZGetP(i*trajGenXYZ.T/n, xyz)      -- stretch out the move requiring less time
    vectorCopy(xyz, 1, world, 1, 3)
    local rslt = err_unrechblpos

    if interpRBT then
      trajGenRBTGetP(i*trajGenRBT.T/n, rbt)    -- stretch out the move requiring less time
      vectorCopy(rbt, 1, theta, 4, 3)
      rslt = inverseWJ(world, theta)
    else -- use Quaternion SLERP
      trajGenQUATGetQ(i*trajGenQUAT.T/n, quat)  -- stretch out the move requiring less time
      vectorCopy(quat, 1, world, 4, 4)
      rslt = inverse(world, theta)
    end

    if rslt~=success then
      print('Error in inverse kinematics: ', rslt)
      break
    end

    -- print out joint command values in radians to the console
    print("[", theta[1], ", ", theta[2], ", ", theta[3], ", ", theta[4], ", ", theta[5], ", ", theta[6], "], ")

  end
end

