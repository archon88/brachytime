# brachytime

A simple Python utility to compute the treatment time (in hours) for an interstitial applicator with wires forming the vertices of equilateral triangles in its midplane, following the Paris system.

A suitable applicator geometry is as shown below, with wire lengths (vertices of the equilateral triangles, denoted by black dots) in millimetres. The basal dose points are denoted by grey dots. OLS quadratic fits were performed to log10(doserate) against log10(distance) to establish a functional relationship between these quantities for wires of length 50 mm, 60 mm, and 70 mm using the data in [AKR.pdf](https://github.com/archon88/brachytime/blob/master/AKR.pdf); linear interpolation is performed for wires of intermediate length using the **interpdose** function. It is assumed that these functional relationships hold for isotopes other than Ir-192. The [WebPlotDigitizer](https://automeris.io/WebPlotDigitizer/) tool was used to convert the data in AKR.pdf to the quadratics in the **calcdoserate** function. The change in source activity over the course of treatment is neglected.

![](diagram.png)

# usage

The treatment time may be computed using expressions such as

```Python
>>> brachytime.calculate_treatment_time(separation=18, upper_lengths=[50, 55, 55], lower_lengths=[70]*4, delivery_akr=0.494, delivery_datetime='2000-11-17T12:00', mid_treatment_datetime='2000-11-22T18:40', halflife=6379000, prescribed_dose=20, printout=True)
```
where the separation (the side of the triangle) and the lengths of the wires in the upper and lower planes are expressed in millimetres, the source air kerma rate at time of delivery is expressed in micrograys per hour per millimetre metre squared, the souce isotope half-life is expressed in seconds, and the prescribed dose is expressed in grays. The source delivery time and the mid-implant time must be passed as strings formatted according to ISO 8601, *i.e.* YYYY-MM-DDTHH:MM. Without loss of generality, it is assumed that the lower row of wires has at least as many wires as the upper row, and that the leftmost wire is in the lower row (since this is taken to be the origin of the coordinate system, for convenience). All wire lengths must be in the range 50 mm to 70 mm.
