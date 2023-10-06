# phaseTurbulenceFields

a multiphase capable version of the OpenFOAM turbulenceFields function object
developed and tested for OFv2306 (with the reactingTwoPhaseEulerFoam solver).

NOTE : this is really dirty and put together fast and little consideration for code quality was made when implementing this. fell free if anyone wants to improve i would be happy to accept PRs !

## Example

add this to your `controlDict`

```c++
    turbulenceFields1
    {
        type            phaseTurbulenceFields;
        libs            (libphaseTurbulenceFields);
        enabled         true;
        writeControl    writeTime;

        phase           "water"; // the name of the phase that you want to calculate the fields for
        fields
        (
            R k // the Reynolds stress tensor R and the turbulent kinetic energy k 
        );
    }
```

the fields are written as `turbulenceProperties:R.water` in this example.

## Building
source your OF environment and run the command `wmake all` in this directory. the `.so` will be created in your `platforms/` directory.

## TODOs:
- add a way to calculate directly the rsm of the velocity
- cleanup the enums and unify the control flow
- more testing
