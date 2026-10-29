# corona3d

The interaction of a coronal mass ejection (CME) with the surrounding corona depends on quantities that vary in all three dimensions, since the electron density sets the height at which a given radio frequency is emitted, the magnetic field and the density together set the Alfvén speed that a shock must exceed, and the topology of the field decides whether accelerated electrons escape or remain trapped on closed loops. Coronagraph and radio images, however, show these processes only in projection, and a three-dimensional picture that combines a coronal model, the CME, and the radio sources remains poorly quantified for most events. corona3d is a Python package that builds this reconstruction for any event from public data, and it follows the method described in Morosan et al. (2022).

In particular, corona3d renders the magnetohydrodynamic (MHD) model MAS/MAST of Predictive Science Inc. (PSI) in three dimensions with PyVista and the Visualization Toolkit (VTK), as a surface of constant electron density coloured by the Alfvén speed together with open and closed field lines. It fits a graduated cylindrical shell (GCS) of the CME to coronagraph images with PyThea, and it places the white-light images of the Large Angle and Spectrometric Coronagraph (LASCO) C2 on the Solar and Heliospheric Observatory (SOHO) and of the COR1 and COR2 coronagraphs on the Solar Terrestrial Relations Observatory (STEREO) as ratio planes in the plane of the sky of their spacecraft. Moreover, the sources of a radio burst are measured in the images of the Nançay Radioheliograph (NRH) and of the Owens Valley Radio Observatory Long Wavelength Array (OVRO-LWA), as blobs and contours, and placed at their heights in a coronal density model on the side of the plane of the sky that faces the observer. The reconstruction is drawn from several observers, such as Earth, a STEREO spacecraft, the solar north pole, or any chosen direction, and one configuration file drives every command, from the download of the data to the figures and the exported layers.

The CME of 2013 October 25, studied by Morosan et al. (2022), is the benchmark against which corona3d is validated, and the CME and type II burst of 2024 May 9 is a new application of the same steps.

![The reconstruction of the benchmark event of 2013 October 25 drawn by corona3d](corona3d_2013-10-25.png)

## Documentation

| Document | Content |
|---|---|
| [User guide](USER_GUIDE.md) | installation, the working directory, the commands and their options, the layers, the viewpoints, the configuration keys, the Python interface, a new event, and the validation against Morosan et al. (2022) |
| [Methods](METHODS.md) | how every layer is computed, which values it uses, and where each value comes from |
| [Limitations](LIMITATIONS.md) | where the results of corona3d can no longer be trusted, and why |
| [Model audit](MODEL_AUDIT.md) | every method and number of corona3d with its source and its test |

The source code, the changes between versions, and the citation file are kept in the [repository](https://github.com/jovian-explorer/corona3d).
