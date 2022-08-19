# REOART 2D and 3D programs
The purpose of this work was to build two ocean-acoustic programs in 2D and 3D able to predict the propagation of sound in a real ocean environment.

The propagation models used are ocean-acoustic ray tracing parametric models based on Fermat’s least time principle and they are composed by systems of algebraic-diﬀerential equations. These systems of equations that model the acoustic trajectories are established through the Euler-Lagrange’s equations.

To study, predict and simulate the propagation of acoustic rays in a real ocean environment is important to have real data about the bathymetry and climatology.
In order to do that we used three diﬀerent databases that give us that kind of information. The databases are GEBCO, WOD18 and TEOS-10.

The GEBCO database enables characterize the ocean bathymetry almost worldwide. In this work it as used the latest version of this database relative to the Portugal oceanic coast.

The WOD18 is a database that contains a variety of data types relevant to oceanography. To calculate the speed of sound in underwater we used the salinity and temperature values from this database.

The TEOS-10 is a software library that contains numerous oceanographic functions for computing relevant ocean variables, and we used it to calculate the speed of sound along the ray trajectory.

Post this, two programs were developed with the ability to study sound propagation in a real ocean environment. The names of the programs are REOART 2D and REOART 3D (Real Environment Ocean-Acoustic Ray Tracing). These programs use the databases mentioned before to withdraw information regarding the ocean topography and climatology to compute the acoustic ray trajectories in a real ocean environment.

Both programs were submitted to several tests and the results were credible and reliable. Not only do the trajectories seem credible, but also the reﬂections on the surface that at the bottom of the ocean seem to satisfy Snell’s law. Given the uniqueness of calculating the speed of sound along the path of the acoustic ray,
these programs could not be compared with other propagation models because of the diﬀerence between the speed of sound.
