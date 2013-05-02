ChConstructor
=============
Building:
---------
make sure you have 'maven' installed then execute

    mvn compile assembly:single

alternatively we have provided a single command Unix script for this

    ./build.sh

and you can find the excutable jar with bundled dependencies in the
target directory.

Running:
--------

    ./chconstructor --help

This will also try to build it if it can't find the jar.

Hardware:
---------

The wrapper script will allocate 20GB java heap memory, and for the german
dataset at least 16GB real memory is recommended; then it will run in
about 5 minutes. With swap it will more likely need hours.

IDE Integration:
----------------
Because we use maven this project can easily be imported into several IDEs
including IntelliJ and Eclipse, just look for 'maven import'
