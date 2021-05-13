Run main.py using Python 3 (tested with Anaconda Python 3.7 on Windows 10). Numpy, scipy, pandas, and matplotlib are required, and requirements.txt describes the versions of all packages used.

Then, interpret the results by hand:
1. First, look at search_results.csv.
   - Discard phases if their cell parameters had to vary too much to produce a good fit.
       -  A good place to find the acceptable range of cell parameters for a given mineral is https://rruff.info/ima/ (check the "cell parameters" box at the top of the page).
   - Discard phases if they have too much of an element not commonly found on Mars.
2. For the remaining phases, look at their plots or CSVs.
   - Discard phases that did not fit the 9.22 A peak well.
   - Discard phases that have significant peaks in places where the Kilmarie pattern does not.
