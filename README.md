# six-axis-robot-library
Library for performing six axis robot forward and inverse kinematics, and trajectory generation (both linear and rotational).
Library supports smooth motion through the wrist singularity of six axis arms. This is done using linear wrist joint interpolation through the singularity while only controlling the end effector x, y, z position. For more details please refer to "DESCRIPTION OF ILLUSTRATIVE EMBODIMENT" section D.2 of this patent: https://patents.google.com/patent/US4680519A/en
Library is capable of generating motion on any 6 axis robot with a spherical wrist.
Library also contains api's for conversion between rotation matrices, quaternions, and Euler angles.
Library uses SLERP (https://en.wikipedia.org/wiki/Slerp) of quaternions for rotational trajectory generation.
For more details on the kinematic model and derivation please refer to appendix A of the thesis found here: https://escholarship.org/uc/item/78h7x9r7

Library has only been tested with lua 5.4.0.
