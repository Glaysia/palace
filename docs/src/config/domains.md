```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# `config["Domains"]`

```json
"Domains":
{
    "Materials":
    [
        ...
    ],
    "CurrentDipole":
    [
        ...
    ],
    "VolumeCurrent":
    [
        ...
    ],
    "Postprocessing":
    {
        "Energy":
        [
            ...
        ],
        "Probe":
        [
            ...
        ]
    }
}
```

with

`"Materials"` :  Array of material properties objects.

`"CurrentDipole"` :  Array of objects for configuring current dipole source excitations.

`"VolumeCurrent"` :  Array of objects for configuring volume current density source
excitations.

`"Postprocessing"` :  Top-level object for configuring domain postprocessing.

`"Energy"` :  Array of objects for postprocessing domain energies.

`"Probe"` :  Array of objects for postprocessing solution field values evaluated at a probe
location in space.

## `domains["Materials"]`

```json
"Materials":
[
    // Material 1
    {
        "Attributes": [<int array>],
        "Permeability": <float> or [<float array>],
        "PermeabilityImag": <float> or [<float array>],
        "MagneticLossTan": <float> or [<float array>],
        "PermeabilityFreq": {
            "Freq": [<float array>],
            "Real": [<float array>],
            "Imag": [<float array>],
            "LossTan": [<float array>]
        },
        "Permittivity": <float> or [<float array>],
        "PermittivityFreq": {
            "Freq": [<float array>],
            "Real": [<float array>],
            "Imag": [<float array>],
            "LossTan": [<float array>]
        },
        "LossTan": <float> or [<float array>],
        "Conductivity": <float> or [<float array>],
        "LondonDepth": <float>,
        "MaterialAxes": [[<array of float array>]]
    },
    // Material 2, 3, ...
    ...
]
```

with

`"Attributes" [None]` :  Integer array of mesh domain attributes for this material.

`"Permeability" [1.0]` :  Relative permeability for this material. Scalar or vector of 3
coefficients corresponding to each of `"MaterialAxes"`.

`"PermeabilityImag" [0.0]` :  Positive magnetic loss component μ″ for complex relative
permeability μ = μ′ - i μ″. Scalar or vector of 3 coefficients corresponding to each of
`"MaterialAxes"`. Specify only one of `"PermeabilityImag"` or `"MagneticLossTan"`.

`"MagneticLossTan" [0.0]` :  Magnetic loss tangent tanδₘ = μ″ / μ′. Scalar or vector of 3
coefficients corresponding to each of `"MaterialAxes"`. Specify only one of
`"MagneticLossTan"` or `"PermeabilityImag"`.

`"PermeabilityFreq" [None]` :  Frequency-dependent scalar relative permeability table for
Driven simulations. `"Freq"` is in GHz and must be strictly increasing. `"Real"` gives μ′.
Use exactly one of `"Imag"` for positive μ″ in μ = μ′ - i μ″ or `"LossTan"` for tanδₘ.
Values are linearly interpolated and solving outside the table range is an error. Frequency-
dependent material tables are not supported together with `"WavePort"` boundaries.

`"Permittivity" [1.0]` : Relative permittivity for this material. Scalar or vector of 3
coefficients corresponding to each of `"MaterialAxes"`.

`"PermittivityFreq" [None]` :  Frequency-dependent scalar relative permittivity table for
Driven simulations. `"Freq"` is in GHz and must be strictly increasing. `"Real"` gives ε′.
Use exactly one of `"Imag"` for positive ε″ in ε = ε′ - i ε″ or `"LossTan"` for tanδ.
Values are linearly interpolated and solving outside the table range is an error. Frequency-
dependent material tables are not supported together with `"WavePort"` boundaries.

`"LossTan" [0.0]` :  Loss tangent for this material. Scalar or vector of 3 coefficients
corresponding to each of `"MaterialAxes"`.

`"Conductivity" [0.0]` :  Electrical conductivity for this material, S/m. Activates Ohmic
loss model in the material domain. Scalar or vector of 3 coefficients corresponding to each
of `"MaterialAxes"`.

`"LondonDepth" [0.0]` :  London penetration depth for this material, specified in mesh
length units. Activates London equations-based model relating superconducting current and
electromagnetic fields in the material domain.

`"MaterialAxes" [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]` : Axes directions for
specification of anisotropic material properties. Required to be unit length and orthogonal.

## `domains["CurrentDipole"]`

```json
"CurrentDipole":
[
    {
        "Index": <int>,
        "Direction": [<float array>],
        "Moment": <float>,
        "Center": [<float array>]
    },
    ...
]
```

with

`"Index" [None]` :  Index of this current dipole source, used in postprocessing output files.

`"Direction" [None]` :  Direction of the Dirac current source specifying the dipole. Axis aligned directions can be specified using
keywords: `"+X"`, `"-X"`, `"+Y"`, `"-Y"`, `"+Z"`, `"-Z"`. The direction can alternatively be specified as a
normalized array of three values, for example `[0.0, 1.0, 0.0]`.

`"Moment" [None]` :  Current dipole moment magnitude, specified in A·m.

`"Center" [None]` :  Floating point array of length equal to the model spatial dimension
specifying the coordinates of the current dipole center position in mesh length units.

## `domains["VolumeCurrent"]`

```json
"VolumeCurrent":
[
    {
        "Index": <int>,
        "Excitation": <int>,
        "Attributes": [<int array>],
        "Direction": [<float array>],
        "CurrentDensity": <float>
    },
    ...
]
```

with

`"Index" [None]` :  Index of this volume current source, used in excitation bookkeeping.

`"Excitation" [None]` :  Driven excitation index for this source. The source contributes
only to the matching excitation.

`"Attributes" [None]` :  Integer array of mesh domain attributes where this source current
density is applied. These attributes must also have corresponding `"Materials"` entries.

`"Direction" [None]` :  Direction of the impressed current density. Axis aligned directions
can be specified using keywords: `"+X"`, `"-X"`, `"+Y"`, `"-Y"`, `"+Z"`, `"-Z"`. The
direction can alternatively be specified as a normalized array of three values, for example
`[0.0, 1.0, 0.0]`.

`"CurrentDensity" [None]` :  Current density magnitude, specified in A/m^2.

## `domains["Postprocessing"]["Energy"]`

```json
"Postprocessing":
{
    "Energy":
    [
        {
            "Index": <int>,
            "Attributes": [<int array>]
        },
        ...
    ]
}
```

with

`"Index" [None]` :  Index of this energy postprocessing domain, used in postprocessing
output files.

`"Attributes" [None]` :  Integer array of mesh domain attributes for this energy
postprocessing domain.

## `domains["Postprocessing"]["Probe"]`

```json
"Postprocessing":
{
    "Probe":
    [
        {
            "Index": <int>,
            "Center": [<float array>]
        },
        ...
    ]
}
```

with

`"Index" [None]` :  Index of this probe, used in postprocessing output files.

`"Center" [None]` :  Floating point array of length equal to the model spatial dimension
specifying the coordinates of this probe in mesh length units.
