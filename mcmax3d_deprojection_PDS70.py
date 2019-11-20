from mcmax3d_analysis.mcmax3d_deprojection import deproject_coordinates

pos_b=(96.8/1000,-147.9/1000)
pos_c=(-233.7/1000,28.8/1000)

print(deproject_coordinates(b=pos_b,c=pos_c,d=113.4,inc=49.7,PA=150))
