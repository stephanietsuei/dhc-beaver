dhc-beaver
==========

Computes the state derivative and output of the DHC-2 "Beaver" aircraft.
This file translates the Simulink model that can be downloaded from
dutchroll.com into a MATLAB function that takes the current state and
controls as input.

It is possible to use this file to compute multiple derivatives and
outputs at once. For example, let

       t = [t1, t2, t3, t4, t5];
       x = [x(t1), x(t2), x(t3), x(t4), x(t5)];
       uaero = [ua(t1), ua(t2), ua(t3), ua(t4), ua(t5)];
       uprop = [up(t1), up(t2), up(t3), up(t4), up(t5)];
       uwind = [uw(t1), uw(t2), uw(t3), uw(t4), uw(t5)];

This function will then compute the following derivative vector in
addition to many outputs.

      dx/dt = [f(x(t1)), f(x(t2)), f(x(t3)), f(x(t4)), f(x(t5))]


Please read the FDC user manual, found at

	https://sourceforge.net/projects/dutchroll/

For more information about the actual flight dynamics.



Licensing
---------

This code is a derivative work of the Flight Dynamics Toolbox, which was
distributed under version 2.1 of the Open Software License, included in the file
license.txt.

Copyright 2013-2014 California Institute of Technology
