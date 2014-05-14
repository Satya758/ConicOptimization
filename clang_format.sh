#!/bin/bash

shopt -s globstar
shopt -s nullglob
shopt -s nocaseglob


for c_file in $(find ./ -name '*.cpp' -or -name '*.hpp')
do
  clang-format-3.4 -style=Google -i "$c_file"
  echo "$c_file"	
done
