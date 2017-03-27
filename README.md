# morse_smale
Python tools for Morse Smale Complex analysis and visualization.

It looks at the critical points of a manifold and partitions the manifold based
on these extrema:

![Example Screenshot](/images/screenshot.png)

# Usage

The common usage will be to create partitions for given sample points and scalar
sample values:

```
import morse_smale

# k_neighbors is the main parameter that may need to be tweaked
result = morse_smale.nn_merged_partition(sample_points, sample_values,
                                         k_neighbors=15)
```

# References

This work is a port of the [msr R package](https://github.com/cran/msr) to
Python.
