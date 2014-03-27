MCNF
====

A program that uses column generation technique to optimize the commodity distribution problem


1. Import the project into Eclipse workspace

2. add the libraries

Expand the MCNF project in the Project Explorer (left hand window) to find the file
“LP.java” from the ’src’ directory. right-click on it, and choose Run As→Run configurations.. 
from the menu that opens. Choose the tab “Environment”, click “New...” and enter
Name: PATH
Value: ${workspace_loc:/MCNF}/jni/${system:OS}_${system:ARCH}
Click “OK”, then “Apply” and then “Run”.
If everything works this should setup and solve a small LP test problem. In particular you
should lookout for the line
lp solve: ifail=0, f=-3.0, LP iters:1
If you see this you are setup for the assignment.

3. Run it. Mcnf.java and Mcnf3.java work properly, Mcnf2 belongs to someone else and I found it in internet.

Copyright
========================================================================
You are free to use it as you like
Send me an email with your comments, good or bad

Christos Giannoukos 2014
chgian@yahoo.com
