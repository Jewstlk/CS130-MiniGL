# Overview
In this project, I implementing a simplified 3D rendering pipeline (with flat shading) in minigl.cpp. This consisted of several parts, introduced through a series of tests:

Vertex and viewing transformations
Rasterization and interpolation
Clipping
Using a z-buffer for hidden surfaces

# Software Framebuffer
Project has a software framebuffer that stores both a 2D array of pixel colors and a 2D array of z (depth) values. My rasterization routine writes into this software framebuffer.

# Rasterization and z-buffer depth test
Implemented a routine that draws a filled triangle specified by three vertices into my software framebuffer. The routine should write its result to the software framebuffer, doing the appropriate z-buffer check.

# Vertex and viewing transformations
Similar to OpenGL, maintains both a projection matrix stack and a modelview stack.

When the user draws a vertex you will apply the modelview and projection ( projection * (modelview * vertex) ) to obtain the transformed geometry. This is followed by a divide by w, and a viewport transform. Stores this transformed geometry for rasterization. The reason for this is that the transformations are part of the current state and will not necessarily persist.

# Tests
This project passes tests 0-25 and can be confirmed by running "./grading-script.sh ."
