#!/usr/bin/env python
"""Wrapper for running panaroo recombination removal from aligned pan-genome"""

# The maintained implementation lives in remove_recombination/.
from remove_recombination.recombination_removal import main

if __name__ == '__main__':
    main()
