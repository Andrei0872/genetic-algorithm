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

def quadratic_fn(a, b, c, x):
  return a * pow(x, 2) + b * x + c

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

def get_chromosome_value(ch: List[int], lb, ub):
  decimal_value = int("".join(list(map(str, ch))), 2)

  return ((ub - lb) * decimal_value) / (pow(2, len(ch)) - 1) + lb

def print_population(config, population: List[List[int]], out_file: TextIOWrapper):
  func = config.function
  
  for i in range(0, len(population)):
    chromosome = population[i]
    ch_value = get_chromosome_value(chromosome, func.lb, func.ub)
    func_value = quadratic_fn(func.a, func.b, func.c, ch_value)

    out_file.write("\t")
    out_file.write("chromosome {}: x = {}; f = {}".format(i + 1, ch_value, func_value))
    out_file.write("\n")

def print_chromosomes_probability(config, population: List[List[int]], out_file: TextIOWrapper):
  func = config.function
  a, b, c = func.a, func.b, func.c
  lb, ub = func.lb, func.ub

  chromosomes_func_values = list(map(lambda ch: quadratic_fn(a, b, c, get_chromosome_value(ch, lb, ub)), population))
  total_func_values = reduce(lambda acc, crt: acc + crt, chromosomes_func_values)

  for i in range(0, len(population)):
    out_file.write("\t")
    out_file.write("chromosome {}: probability of {}".format(i + 1, chromosomes_func_values[i] / total_func_values))
    out_file.write("\n")

def print_selection_intervals(config, population: List[List[int]], out_file: TextIOWrapper):
  func = config.function
  a, b, c = func.a, func.b, func.c
  lb, ub = func.lb, func.ub

  chromosomes_func_values = list(map(lambda ch: quadratic_fn(a, b, c, get_chromosome_value(ch, lb, ub)), population))
  total_func_values = reduce(lambda acc, crt: acc + crt, chromosomes_func_values)
  probabilities = list(map(lambda ch: ch / total_func_values, chromosomes_func_values))
  selection_intervals_points = [0] + list(accumulate(probabilities, lambda x, y: x + y))

  for i in range(0, len(selection_intervals_points) - 1):
    out_file.write("\t")
    out_file.write("interval {}: [{}, {})".format(i + 1, selection_intervals_points[i], selection_intervals_points[i + 1]))
    out_file.write("\n")


if __name__ == "__main__":
  out_file = open(OUTPUT_FILE_NAME, "w")

  config = read_config(INPUT_FILE_NAME)

  print(config)
  
  chromosome_length = get_chromosome_length(config.function.lb, config.function.ub, config.precision)
  initial_pop = generate_population(config.nr_individuals, chromosome_length)

  out_file.write("1. Initial population: \n")
  print_population(config, initial_pop, out_file)
  out_file.write(SECTION_SEPARATOR)

  out_file.write("2. Probability of selection for each chromosome: \n")
  print_chromosomes_probability(config, initial_pop, out_file)
  out_file.write(SECTION_SEPARATOR)

  out_file.write("3. Selection intervals: \n")
  print_selection_intervals(config, initial_pop, out_file)
  out_file.write(SECTION_SEPARATOR)
