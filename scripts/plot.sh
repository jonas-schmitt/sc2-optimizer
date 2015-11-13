#!/bin/bash
gnuplot ./scripts/plot_avg_damage.gnuplot
gnuplot ./scripts/plot_avg_health.gnuplot
gnuplot ./scripts/plot_stdev_damage.gnuplot
gnuplot ./scripts/plot_stdev_health.gnuplot
gnuplot ./scripts/plot_paths.gnuplot
inkscape --export-area-drawing --export-png=avg_damage.png --export-dpi=300 avg_damage.svg
inkscape --export-area-drawing --export-png=avg_health.png --export-dpi=300 avg_health.svg
inkscape --export-area-drawing --export-png=stdev_damage.png --export-dpi=300 stdev_damage.svg
inkscape --export-area-drawing --export-png=stdev_health.png --export-dpi=300 stdev_health.svg
#inkscape --export-area-drawing --export-png=paths.png --export-dpi=300 paths.svg
