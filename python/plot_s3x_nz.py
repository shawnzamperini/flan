"""
The old name of plot_multizone.py, which now plots any multi-zone run
(SOLEDGE3X or HERMES) and any output or background variable. This runs it
with the same arguments, and keeps

	from plot_s3x_nz import plot_s3x_nz, movie_s3x_nz

working.
"""

from plot_multizone import (load_zones, make_norm, draw_zones,  # noqa: F401
	plot_multizone, movie_multizone, plot_s3x_nz, movie_s3x_nz, main)

if __name__ == "__main__":
	main()
