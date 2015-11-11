#!/bin/bash
gnuplot plot_processors.gnuplot
gnuplot plot_nodes.gnuplot
inkscape --export-area-drawing --export-png=weak_scaling_processors.png --export-dpi=300 weak_scaling_processors.svg
inkscape --export-area-drawing --export-png=weak_scaling_nodes.png --export-dpi=300 weak_scaling_nodes.svg
