*Not really tested for MAC implementation.

*For Linux, path to elastix binaries shuold be added to bashrc

*if you get an error that the libraries cannot be found, you may need to export them explicitly - or add an additional line to bashrc

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/PATH_TO_SEEVR/seeVR-main/seeVR-main/registration/elastix/linux/lib

ALSO - be sure that elastix and transformix can be executed as a program - right click properties and check in the permissions tab