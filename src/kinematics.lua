
-------------------------------------------------------------------------------------------------------------
-- Redefine math library 

pi    = math.pi
rad   = math.rad
deg   = math.deg
sin   = math.sin
cos   = math.cos
tan   = math.tan
asin  = math.asin
acos  = math.acos
atan  = math.atan
atan2 = math.atan2
sqrt  = math.sqrt
abs   = math.abs
min   = math.min
max   = math.max

-------------------------------------------------------------------------------------------------------------
-- Robot model (MH5 parameters)
a0 = 88
a1 = 310
a2 = 40
d3 = -305
a4 = 0
d4 = 0
d5 = -80

minTheta = {rad(-170), rad(-65), rad(-70), rad(-180), rad(-135), rad(-180), rad(-136)}
maxTheta = {rad( 170), rad(150), rad(190), rad( 180), rad( 135), rad( 180), rad( 255)}

-------------------------------------------------------------------------------------------------------------
-- Robot model (MHJ parameters)
-- a0= 0
-- a1= 275
-- a2= 0
-- d3= -270
-- a4= 0
-- d4= 0
-- d5= -63

-- minTheta = {rad(-160), rad(-90), rad(-45), rad(-180), rad(-130), rad(-180), rad(-105)}   -- DATASHEET HAS TYPO!!! -290/+105 should be -105/+290
-- maxTheta = {rad( 160), rad(110), rad(210), rad( 180), rad( 130), rad( 180), rad( 290)}   -- DATASHEET HAS TYPO!!! -290/+105 should be -105/+290

-------------------------------------------------------------------------------------------------------------
-- Error codes and tolerances

success = 0
err_unrechblpos = 1
err_unrechblorn = 2
err_joint_limit = 3
err_singularity = 4
err_bifurcation = 5

-- bifurcation config flags
bif1 = 1 -- front/back reach
bif2 = 2 -- upper/lower arm
bif3 = 4 -- Flip/No Flip wrist

Encoder_Resolution = 16777216  -- assume 24 bit encoder
Typical_LinkLength = 1000

tol_ang = 2*pi/Encoder_Resolution
tol_lin = tol_ang*Typical_LinkLength
tol_ang_sng = 1e-15
tol_lin_sng = tol_ang_sng*Typical_LinkLength


theta_offset = {0, -pi/2, 0, 0, 0, 0} -- offset between DH parameters and zero home position

-------------------------------------------------------------------------------------------------------------
-- Pre-compute constant terms

d4d4 = d4*d4
d3d3 = d3*d3
a0a4 = a0+a4
a2a2 = a2*a2
kk3 = a0a4*a0a4-a1*a1-a2a2-d3d3-d4d4
theta2_max_reach = atan(a2/d3)+pi/2
tol_lin_sng2 = tol_lin_sng*tol_lin_sng

-------------------------------------------------------------------------------------------------------------

tol_gimbal = 1e-14

-------------------------------------------------------------------------------------------------------------
-- Global storage for vectors and matrices

theta  = {0, 0, 0, 0, 0, 0}
theta1 = {0, 0, 0, 0, 0, 0}
theta2 = {0, 0, 0, 0, 0, 0}

world  = {0, 0, 0, 0, 0, 0, 0}
world1 = {0, 0, 0, 0, 0, 0, 0}
world2 = {0, 0, 0, 0, 0, 0, 0}

xyz  = {0, 0, 0}
xyz1 = {0, 0, 0}
xyz2 = {0, 0, 0}

rbt  = {0, 0, 0}
rbt1 = {0, 0, 0}
rbt2 = {0, 0, 0}

quat  = {0, 0, 0, 0}
quat1 = {0, 0, 0, 0}
quat2 = {0, 0, 0, 0}

abg  = {0, 0, 0}

R = { { 0, 0, 0 },
      { 0, 0, 0 },
      { 0, 0, 0 } }

A = { { 0, 0, 0, 0 },
      { 0, 0, 0, 0 },
      { 0, 0, 0, 0 },
      { 0, 0, 0, 1 } }

-------------------------------------------------------------------------------------------------------------

function vectorCopy(src, srcIdx, dst, dstIdx, size)
  for i=0,size-1 do
    dst[dstIdx + i] = src[srcIdx+i]
  end
end

function matrixCopy(src, dst, rows, cols)
  for i=1,rows do
    for j=1,cols do
       dst[i][j] = src[i][j]
     end
   end
end

-------------------------------------------------------------------------------------------------------------

function constrain(theta, minimum, maximum)
  if theta<minimum then
    theta = theta+2*pi
  end
  if theta>maximum then
    theta = theta-2*pi
  end
  return theta
end

function normalizeQuat(q)
  local c = 1/sqrt(q[1]*q[1] +q[2]*q[2] +q[3]*q[3] +q[4]*q[4])
  q[1] = q[1]*c
  q[2] = q[2]*c
  q[3] = q[3]*c
  q[4] = q[4]*c
end

-------------------------------------------------------------------------------------------------------------

function eulerToR(abg, R)
  local alpha = abg[1]
  local betta = abg[2]
  local gamma = abg[3]

  local salpha = sin(alpha)
  local calpha = cos(alpha)

  local sbetta = sin(betta)
  local cbetta = cos(betta)

  local sgamma = sin(gamma)
  local cgamma = cos(gamma)

  R[1][1] =  cbetta*cgamma   R[1][2] = -calpha*sgamma+salpha*sbetta*cgamma   R[1][3] =  salpha*sgamma+calpha*sbetta*cgamma
  R[2][1] =  cbetta*sgamma   R[2][2] =  calpha*cgamma+salpha*sbetta*sgamma   R[2][3] = -salpha*cgamma+calpha*sbetta*sgamma
  R[3][1] = -sbetta          R[3][2] =                salpha*cbetta          R[3][3] =                calpha*cbetta
end

function RToEuler(R, abg)
  local alpha
  local betta
  local gamma

  if abs(abs(R[3][1])-1)>tol_gimbal then
    local betta = asin(-R[3][1])
    local cbetta = cos(betta)
    local alpha = atan2(R[3][2]/cbetta, R[3][3]/cbetta)
    local gamma = atan2(R[2][1]/cbetta, R[1][1]/cbetta)
  else
    gamma = 0
    if R[3][1]<0 then
      betta = pi/2
      alpha = gamma+atan2(R[1][2], R[1][3])
    else
      betta = -pi/2
      alpha = -gamma+atan2(-R[1][2], -R[1][3])
    end
  end
  
  abg[1] = alpha
  abg[2] = betta
  abg[3] = gamma
end

-------------------------------------------------------------------------------------------------------------

function eulerToQuat(abg, q)
  local alpha = abg[1]
  local betta = abg[2]
  local gamma = abg[3]

  local sa = sin(alpha/2)
  local sb = sin(betta/2)
  local sg = sin(gamma/2)
  local ca = cos(alpha/2)
  local cb = cos(betta/2)
  local cg = cos(gamma/2)

  q[1] = ca*cb*cg+sa*sb*sg
  q[2] = sa*cb*cg-ca*sb*sg
  q[3] = ca*sb*cg+sa*cb*sg
  q[4] = ca*cb*sg-sa*sb*cg
end

function quatToEuler(q, abg)
  local alpha
  local betta
  local gamma

  local ysqr = q[3]*q[3]

  local t0 = 2.0*(q[1]*q[2]+q[3]*q[4])
  local t1 = 1.0-2.0*(q[2]*q[2]+ysqr)
  local t2 = 2.0*(q[1]*q[3]-q[4]*q[2])
  

  if t2 > 1.0 then
    t2 = 1.0
  end

  if t2 < -1.0 then
    t2 = -1.0
  end

  local t3 = 2.0*(q[1]*q[4]+q[2]*q[3])
  local t4 = 1.0-2.0*(ysqr+q[4]*q[4])

  if abs(abs(t2)-1)>tol_gimbal then
    alpha = atan2(t0, t1)
    betta = asin(t2)
    gamma = atan2(t3, t4)
  else
    gamma = 0
    local t5 = 2*(q[2]*q[3]-q[4]*q[1])
    local t6 = 2*(q[2]*q[4]+q[3]*q[1])
    if t2>0 then
      betta = pi/2
      alpha = gamma+atan2(t5, t6)
    else
      betta = -pi/2
      alpha = -gamma+atan2(-t5, -t6)
    end
  end
  
  abg[1] = alpha
  abg[2] = betta
  abg[3] = gamma
end

-------------------------------------------------------------------------------------------------------------

function quatToR(q, R)

  local qw = q[1]
  local qx = q[2]
  local qy = q[3]
  local qz = q[4]

  local qx2 = qx*qx
  local qy2 = qy*qy
  local qz2 = qz*qz

  local qwqx = qw*qx
  local qwqy = qw*qy
  local qwqz = qw*qz

  local qxqy = qx*qy
  local qxqz = qx*qz

  local qyqz = qy*qz

  R[1][1] =  1 -2*qy2 -2*qz2   R[1][2] = 2*(qxqy - qwqz)   R[1][3] =  2*(qxqz + qwqy)
  R[2][1] =  2*(qxqy + qwqz)   R[2][2] = 1 -2*qx2 -2*qz2   R[2][3] =  2*(qyqz - qwqx)
  R[3][1] =  2*(qxqz - qwqy)   R[3][2] = 2*(qyqz + qwqx)   R[3][3] =  1 -2*qx2 -2*qy2
end

function RToQuat(R, q)
  local qw
  local qx
  local qy
  local qz
  local S

  local tr = R[1][1]+R[2][2]+R[3][3]
  if tr>0 then
    S = sqrt(tr+1.0)*2
    qw = 0.25*S
    qx = (R[3][2]-R[2][3])/S
    qy = (R[1][3]-R[3][1])/S
    qz = (R[2][1]-R[1][2])/S
  elseif R[1][1]>R[2][2] and R[1][1]>R[3][3] then
    S = sqrt(1.0+R[1][1]-R[2][2]-R[3][3])*2
    qw = (R[3][2]-R[2][3])/S
    qx = 0.25*S
    qy = (R[1][2]+R[2][1])/S
    qz = (R[1][3]+R[3][1])/S
  elseif R[2][2]>R[3][3] then
    S = sqrt(1.0+R[2][2]-R[1][1]-R[3][3])*2
    qw = (R[1][3]-R[3][1])/S
    qx = (R[1][2]+R[2][1])/S
    qy = 0.25*S
    qz = (R[2][3]+R[3][2])/S
  else
    S = sqrt(1.0+R[3][3]-R[1][1]-R[2][2])*2
    qw = (R[2][1]-R[1][2])/S
    qx = (R[1][3]+R[3][1])/S
    qy = (R[2][3]+R[3][2])/S
    qz = 0.25*S
  end
  
  q[1] = qw
  q[2] = qx
  q[3] = qy
  q[4] = qz
end

-------------------------------------------------------------------------------------------------------------

function worldToMatrix(world, A)
  vectorCopy(world, 4, quat, 1, 4)
  quatToR(quat, R)
  matrixCopy(R, A, 3, 3)

  A[1][4] = world[1]
  A[2][4] = world[2]
  A[3][4] = world[3]
end

function matrixToWorld(A, world)
  world[1] = A[1][4]
  world[2] = A[2][4]
  world[3] = A[3][4]

  matrixCopy(A, R, 3, 3)
  RToQuat(R, quat)
  vectorCopy(quat, 1, world, 4, 4)
end

-------------------------------------------------------------------------------------------------------------

function constrainthetaPos(theta)
  for x=1,6 do
    if x==3 then
      theta[x] = constrain(theta[x],  -pi/2, 3*pi/2)
    else
      theta[x] = constrain(theta[x], -pi, pi)
    end
  end
end

function verifythetaLimits(theta)
  for x=1, 6 do
    if theta[x]<minTheta[x]-tol_ang or theta[x]>maxTheta[x]+tol_ang then
      return err_joint_limit
    end
  end

  if theta[3]-theta[2]<minTheta[7]-tol_ang or theta[3]-theta[2]>maxTheta[7]+tol_ang then
    return err_joint_limit
  end
  
  return success
end

-------------------------------------------------------------------------------------------------------------

function ModelToDH(theta)   -- convert from kinematic model joint position space to DH model joint position space
  for x=1,6 do
    theta[x] = theta[x]+theta_offset[x]
  end
end

function DHToModel(theta)
  for x=1,6 do
    theta[x] = theta[x]-theta_offset[x]
  end
end

-------------------------------------------------------------------------------------------------------------

function thetaToConfig(theta)
  local config = 0

  ModelToDH(theta)
  
  if a0+a4+a1*cos(theta[2])+a2*cos(theta[3]-theta[2])-d3*sin(theta[3]-theta[2])<0 then
    config = config+bif1
  end
  
  if theta[3]>theta2_max_reach then
    config = config+bif2
  end

  if theta[5]>0 then
    config = config+bif3
  end

  DHToModel(theta)
  
  return config
end

function chooseCloserWristConfig(theta1, theta2)

  local dR1 = theta2[4] - theta1[4]
  local dB1 = theta2[5] - theta1[5]
  local dT1 = theta2[6] - theta1[6]

  local dR2 = constrain(theta2[4] + 2*pi, -2*pi, 2*pi) - theta1[4]
  local dB2 =          -theta2[5]                      - theta1[5]
  local dT2 = constrain(theta2[6] + 2*pi, -2*pi, 2*pi) - theta1[6]

  local dis1 = 7*dR1*dR1 + 3*dB1*dB1 + 1*dT1*dT1  -- use larger weights for rotation distances of R, then B, then T joints, due to their inertias being as such
  local dis2 = 7*dR2*dR2 + 3*dB2*dB2 + 1*dT2*dT2

  if dis2 < dis1 then
    theta2[4] = constrain(theta2[4] + 2*pi, -2*pi, 2*pi)
    theta2[5] =          -theta2[5]
    theta2[6] = constrain(theta2[6] + 2*pi, -2*pi, 2*pi)
  end
end

function constainRJointInWristSing(theta1, theta2)
  if abs(theta2[5]) < tol_ang_sng then
    theta2[6] = theta2[6] + (theta2[4] - theta1[4])
    theta2[4] = theta2[4] - (theta2[4] - theta1[4]) -- ie. theta2[4] = theta1[4]

    theta2[4] = constrain(theta2[4], -2*pi, 2*pi)
    theta2[6] = constrain(theta2[6], -2*pi, 2*pi)
  end
end

function isWristSingInPath(theta1, theta2)
  return  abs(theta1[5]) < pi/6 or abs(theta2[5]) < pi/6 or theta1[5]*theta2[5] < 0
end

function solveSinCosTheta(k1, k2, k3, theta, i)

  if abs(k1)<tol_lin_sng and abs(k2)<tol_lin_sng then
    return err_singularity
  end

  local somega = k3/sqrt(k1*k1 + k2*k2)

  if somega > 1 then
    return err_unrechblpos
  end

  local phi = atan2(-k2, k1)
  local omega = asin(somega)

  local sol1 = phi + omega
  local sol2 = phi + 2*pi - omega

  sol1 = constrain(sol1, -2*pi, 2*pi)
  sol2 = constrain(sol2, -2*pi, 2*pi)

  -- Always pick the solution closer to the last position
  if abs(theta[i] - sol1) < abs(theta[i] - sol2) then
    theta[i] = sol1
  else
    theta[i] = sol2
  end

  return success
end

function forward(theta, world)
  local rslt = verifythetaLimits(theta)

  ModelToDH(theta)

  local theta2_1 = theta[3]-theta[2]
  local stheta2_1 = sin(theta2_1)
  local ctheta2_1 = cos(theta2_1)
  local ctheta0 = cos(theta[1])
  local stheta0 = sin(theta[1])
  local ctheta1 = cos(theta[2])
  local stheta1 = sin(theta[2])
  local ctheta3 = cos(theta[4])
  local stheta3 = sin(theta[4])
  local ctheta4 = cos(theta[5])
  local stheta4 = sin(theta[5])
  local ctheta5 = cos(theta[6])
  local stheta5 = sin(theta[6])

  local ctheta0ctheta2_1 = ctheta0*ctheta2_1
  local stheta0ctheta2_1 = stheta0*ctheta2_1
  local stheta0stheta2_1 = stheta0*stheta2_1
  local ctheta0stheta2_1 = ctheta0*stheta2_1
  local stheta4ctheta5 = stheta4*ctheta5
  local stheta3stheta4 = stheta3*stheta4
  local ctheta3ctheta4 = ctheta3*ctheta4
  local stheta4stheta5 = stheta4*stheta5
  local ctheta3stheta4 = ctheta3*stheta4
  local d5ctheta3stheta4 = d5*ctheta3stheta4
  local d5stheta3stheta4 = d5*stheta3stheta4

  local x = ctheta3ctheta4*ctheta5-stheta3*stheta5
  local y = stheta3*ctheta4*ctheta5+ctheta3*stheta5
  local z = ctheta3ctheta4*stheta5+stheta3*ctheta5
  local w = stheta3*ctheta4*stheta5-ctheta3*ctheta5
  local m = d5*ctheta4+d3
  local n = a0+a1*ctheta1+a2*ctheta2_1

  local ux = ctheta0ctheta2_1*x-stheta0*y-ctheta0stheta2_1*stheta4ctheta5
  local uy = stheta0ctheta2_1*x+ctheta0*y-stheta0stheta2_1*stheta4ctheta5
  local uz = stheta2_1*x+ctheta2_1*stheta4ctheta5
  local vx = ctheta0ctheta2_1*z-stheta0*w-ctheta0stheta2_1*stheta4stheta5
  local vy = stheta0ctheta2_1*z+ctheta0*w-stheta0stheta2_1*stheta4stheta5
  local vz = stheta2_1*z+ctheta2_1*stheta4stheta5
  local wx = ctheta0ctheta2_1*ctheta3stheta4-stheta0*stheta3stheta4+ctheta0stheta2_1*ctheta4
  local wy = stheta0ctheta2_1*ctheta3stheta4+ctheta0*stheta3stheta4+stheta0stheta2_1*ctheta4
  local wz = stheta2_1*ctheta3stheta4-ctheta2_1*ctheta4
  local qx = -ctheta0ctheta2_1*d5ctheta3stheta4+stheta0*d5stheta3stheta4-ctheta0stheta2_1*m+ctheta0*n
  local qy = -stheta0ctheta2_1*d5ctheta3stheta4-ctheta0*d5stheta3stheta4-stheta0stheta2_1*m+stheta0*n
  local qz = -stheta2_1*d5ctheta3stheta4+ctheta2_1*m-a1*stheta1+a2*stheta2_1

  A[1][1] = ux   A[1][2] = vx   A[1][3] = wx   A[1][4] = qx
  A[2][1] = uy   A[2][2] = vy   A[2][3] = wy   A[2][4] = qy
  A[3][1] = uz   A[3][2] = vz   A[3][3] = wz   A[3][4] = qz

  DHToModel(theta)

  matrixToWorld(A, world)
  return rslt
end

function inverse(world, theta)
  ModelToDH(theta)

  local thetaPrev = {theta[1], theta[2], theta[3], theta[4], theta[5], theta[6]}

  worldToMatrix(world, A)
  
  local ux = A[1][1]
  local uy = A[2][1]
  local uz = A[3][1]
  local vx = A[1][2]
  local vy = A[2][2]
  local vz = A[3][2]
  local wx = A[1][3]
  local wy = A[2][3]
  local wz = A[3][3]
  local qx = A[1][4]
  local qy = A[2][4]
  local qz = A[3][4]

  local px = qx+d5*wx
  local py = qy+d5*wy
  local pz = qz+d5*wz

  local k1 = px
  local k2 = -py
  local k3 = d4

  --bif1, choose closer configuration
  if solveSinCosTheta(k1, k2, k3, theta, 1) ~= success then
    return err_unrechblpos
  end

  local stheta0 = sin(theta[1])
  local ctheta0 = cos(theta[1])

  k1 = -2*a1*d3
  k2 = 2*a1*a2
  k3 = px*px+py*py+pz*pz-2*px*a0a4*ctheta0-2*py*a0a4*stheta0+kk3
  
  --bif2, choose closer configuration
  if solveSinCosTheta(k1, k2, k3, theta, 3) ~= success then
    return err_unrechblpos
  end
  
  local stheta2 = sin(theta[3])
  local ctheta2 = cos(theta[3])
  
  local u1 = a1+a2*ctheta2-d3*stheta2
  local v1 = a2*stheta2+d3*ctheta2
  local y1 = px*ctheta0+py*stheta0-a0a4
  local u2 = -a2*stheta2-d3*ctheta2
  local v2 = a1+a2*ctheta2-d3*stheta2
  local y2 = -pz
  
  if abs(v2*u1-v1*u2)<tol_lin_sng2 then
    return err_singularity
  end
  theta[2] = atan2(u1*y2-u2*y1, v2*y1-v1*y2)

  local theta2_1 = theta[3]-theta[2]
  local stheta2_1 = sin(theta2_1)
  local ctheta2_1 = cos(theta2_1)
  

  theta[5] = acos(wx*ctheta0*stheta2_1+wy*stheta0*stheta2_1-wz*ctheta2_1)

  if abs(theta[5])>tol_ang_sng then
    theta[4] = atan2(-wx*stheta0+wy*ctheta0, wx*ctheta0*ctheta2_1+wy*stheta0*ctheta2_1+wz*stheta2_1)
    theta[6] = atan2(-vx*ctheta0*stheta2_1-vy*stheta0*stheta2_1+vz*ctheta2_1, -ux*ctheta0*stheta2_1-uy*stheta0*stheta2_1+uz*ctheta2_1)
    
    --bif3, choose closer configuration
    chooseCloserWristConfig(thetaPrev, theta)
  else -- we are in wrist singularity, so keep theta[4] unchanged
    theta[6] = atan2(-ux*stheta0 + uy*ctheta0, ux*ctheta0*ctheta2_1 + uy*stheta0*ctheta2_1 + uz*stheta2_1) - theta[4]
  end

  DHToModel(theta)
  constrainthetaPos(theta)
  return verifythetaLimits(theta)
end

function inverseWJ(world, theta)

  ModelToDH(theta)

  local qx = world[1]
  local qy = world[2]
  local qz = world[3]

  local ctheta3 = cos(theta[4])
  local stheta3 = sin(theta[4])
  local ctheta4 = cos(theta[5])
  local stheta4 = sin(theta[5])

  local m = ctheta3*stheta4
  local n = stheta3*stheta4
  local l = d5*ctheta4+d3

  local k1 = qx
  local k2 = -qy
  local k3 = n*d5

  if solveSinCosTheta(k1, k2, k3, theta, 1) ~= success then
    return err_unrechblpos
  end

  local ctheta0 = cos(theta[1])
  local stheta0 = sin(theta[1])

  local r = ctheta0*qx + stheta0*qy - a0

  k1 = 2*qz*a1
  k2 = -2*r*a1
  k3 = (a2 - m*d5)*(a2 - m*d5) + l*l - r*r-qz*qz - a1*a1

  if solveSinCosTheta(k1, k2, k3, theta, 2) ~= success then
    return err_unrechblpos
  end

  local ctheta1 = cos(theta[2])
  local stheta1 = sin(theta[2])

  local u1 = a2-m*d5
  local v1 = -l
  local y1 = r-a1*ctheta1
  local u2 =l
  local v2 = a2-m*d5
  local y2 = qz+a1*stheta1
  
  if abs(v2*u1-v1*u2)<tol_lin_sng2 then
    return err_singularity
  end

  local theta2_1 = atan2(u1*y2-u2*y1, v2*y1-v1*y2)
  
  theta[3] = theta2_1 + theta[2]

  DHToModel(theta)
  constrainthetaPos(theta)
  return verifythetaLimits(theta)
end

