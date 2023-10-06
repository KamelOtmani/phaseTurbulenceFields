# phaseTurbulenceFields

a multiphase capable version of the OpenFOAM turbulenceFields function object
developed and tested for OFv2306

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