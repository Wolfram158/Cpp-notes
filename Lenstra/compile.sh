#!/bin/bash

g++ -std=c++17 -I . -pthread main.cpp -o factor -lgmpxx -lgmp -pedantic -Wall -Werror
