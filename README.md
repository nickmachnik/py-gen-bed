# py-gen-bed

Generate a random .bed file:

```python
from gen_bed_py import rand_bed_file

bed_output_file = "./my.bed"
number_of_individuals = 1000
number_of_markeres = 10000

rand_bed_file(number_of_individuals, number_of_markers, bed_output_file)
```

## DISCLAIMER

None of this code is tested, use at your own risk in production.