awk '{if (NR%2==0) printf("{ frequency %10.5e values { kw = %10.5e 1.0e+0 log 1  } }\n", $1, $2)}'
