# morse_smale
Python tools for Morse Smale Complex analysis and visualization.

It looks at the critical points of a manifold and partitions the manifold based
on these extrema.

![Example Screenshot](/images/screenshot.png)

# Usage

The common usage will be to create partitions for given sample points and scalar
sample values.

```
import morse_smale

# k_neighbors is the main parameter that may need to be tweaked
result = morse_smale.nn_merged_partition(sample_points, sample_values,
                                         k_neighbors=15)
```

See `examples/` for more.

# References

This work is a port of the [msr R package](https://github.com/cran/msr) to
Python. It builds on an interface to the C++ Morse-Smale complex computation
first described in:

[Gerber S, Bremer PT, Pascucci V, Whitaker R (2010). "Visual Exploration of High 
Dimensional Scalar Functions." IEEE Transactions on Visualization and Computer 
Graphics, 16(6), 1271–1280.](https://www.ncbi.nlm.nih.gov/pubmed/20975167)
