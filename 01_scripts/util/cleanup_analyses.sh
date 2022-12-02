#!/bin/bash
# Remove any result and temporary file from folders 05_trimmed and up

# Cleanup temporary data and result directories
rm 05_trimmed/* 2>/dev/null
rm 06_extracted/* 2>/dev/null
rm 07_aligned/* 2>/dev/null
