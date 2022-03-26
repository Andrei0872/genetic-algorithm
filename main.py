from functools import reduce
from io import TextIOWrapper
from itertools import accumulate
from math import ceil, log2
from types import SimpleNamespace
import random
from typing import List

INPUT_FILE_NAME = "input.txt"
OUTPUT_FILE_NAME = "output.txt"
SECTION_SEPARATOR = "\n" + 50 * "=" + "\n\n"

def read_config(input_file_name) -> SimpleNamespace:
  input_file = open(input_file_name, "r")
  
  config = SimpleNamespace()

  config.nr_individuals = int(input_file.readline())

  config.function = SimpleNamespace()
  config.function.lb = int(input_file.readline())
  config.function.ub = int(input_file.readline())
  config.function.a = int(input_file.readline())
  config.function.b = int(input_file.readline())
  config.function.c = int(input_file.readline())

  config.precision = int(input_file.readline())
  config.prob_crossover = float(input_file.readline())
  config.prob_mutation = float(input_file.readline())

  config.nr_phases = int(input_file.readline())

  input_file.close()
  return config

def generate_chromosome(len):
  chromosome = []

  for _ in range(0, len):
    gene = round(random.random())
    chromosome.append(gene)
  
  return chromosome

def generate_population(nr_individuals, chromosome_length):
  population = []

  for i in range(0, nr_individuals):
    chromosome = generate_chromosome(chromosome_length)
    population.append(chromosome)
  
  return population

def get_chromosome_length(lb, ub, precision):
  return ceil(log2((ub - lb) * pow(10, precision)))

if __name__ == "__main__":
  out_file = open(OUTPUT_FILE_NAME, "w")

  config = read_config(INPUT_FILE_NAME)

  print(config)
  
  chromosome_length = get_chromosome_length(config.function.lb, config.function.ub, config.precision)
  initial_pop = generate_population(config.nr_individuals, chromosome_length)

  out_file.write("1. Initial population: \n")
  print_population(config, initial_pop, out_file)
  out_file.write(SECTION_SEPARATOR)

