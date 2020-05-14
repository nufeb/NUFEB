#include "shapes.inc"
// Right-handed coordinate system in which the z-axis points upwards
camera {
location <0, 60e-2, 0>
sky z
right 0.24*x*image_width/image_height
up 0.24*y
look_at <0, 0, 0>
look_at 0
rotate <0, 0, 20>
}
// Create simualtion domain
#declare b_x=0.05;
#declare b_y=0.02;
#declare b_z=0.05;
object{ 
Wire_Box(<b_x, b_y, b_z>, <-b_x, -b_y, -b_z>, 0.0002, 0) 
 texture{ pigment{ color rgb<1,1,1>}}
}
// White background
background{rgb 1}
// Two lights with slightly different colors
light_source{<-10e-2,80e-2,0e-2> color rgb <0.77,0.75,0.75>}
light_source{<10e-2,80e-2,0e-2> color rgb <0.57,0.55,0.55>}
// Radius of the Voronoi cell network
#declare r=0.05;
// Radius of the particles
#declare s=0.6;
// Particles
union{
#include "data_p.pov"
}
