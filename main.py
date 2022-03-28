from copy import deepcopy
from enum import Enum
from functools import reduce
from io import TextIOWrapper
from itertools import accumulate
from math import ceil, log2
from types import SimpleNamespace
import random
from typing import List, Tuple

INPUT_FILE_NAME = "input.txt"
OUTPUT_FILE_NAME = "output.txt"
SECTION_SEPARATOR = "\n" + 50 * "=" + "\n\n"

def print_if_true(discriminant, fn):
  if discriminant:
    fn()

class MutationTypes(Enum):
  RARE = 1
  EACH_GENE = 2

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
    selection_intervals.append((selection_intervals_points[i], selection_intervals_points[i + 1]))
    
    if out_file != None:
      out_file.write("\t")
      out_file.write("interval {}: [{}, {})".format(i + 1, selection_intervals_points[i], selection_intervals_points[i + 1]))
      out_file.write("\n")
  
  return selection_intervals

def find_chromosome_idx_in_interval(intervals: List[Tuple[float, float]], prob_chromosome):
  start, end = 0, len(intervals) - 1

  while start <= end:
    mid = (start + end) // 2

    (lb, ub) = intervals[mid]
    if lb <= prob_chromosome < ub:
      return mid
    
    if prob_chromosome > ub:
      start = mid + 1
    elif prob_chromosome < lb:
      end = mid - 1

def select_from_population(population: List[List[int]], selection_intervals: List[Tuple[float, float]], out_file: TextIOWrapper = None) -> List[List[int]]:
  selected_population = []
  
  for _ in population:
    prob_chromosome = random.uniform(0, 1)
    chromosome_idx = find_chromosome_idx_in_interval(selection_intervals, prob_chromosome)

    if out_file != None:
      out_file.write("\t")
      out_file.write("prob_chromosome = {}; selecting chromosome {}".format(prob_chromosome, chromosome_idx + 1))
      out_file.write("\n")

    selected_population.append(deepcopy(population[chromosome_idx]))
  
  return selected_population

def split_based_on_crossover_prob(selected_population: List[List[int]], prob_crossover, out_file: TextIOWrapper = None) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
  in_crossover = []
  out_crossover = []
  
  for idx in range(0, len(selected_population)):
    prob_chromosome = random.uniform(0, 1)
    chromosome = selected_population[idx]
    
    if out_file != None:
      out_file.write("\t")
      out_file.write("chromosome {}: {}; prob_chromosome = {}".format(idx + 1, "".join(list(map(str, chromosome))), prob_chromosome))

    if prob_chromosome >= prob_crossover:
      out_crossover.append((idx, chromosome))
      if out_file != None:
        out_file.write("\n")

      continue
    
    if out_file != None:
      out_file.write(" < {} ===> selected\n".format(prob_crossover))

    in_crossover.append((idx, chromosome))

  return (in_crossover, out_crossover)

def chromosomes_crossover(ch1: List[int], ch2: List[int], breakpoint) -> Tuple[List[int], List[int]]:
  after_breakpoint_ch1 = ch1[breakpoint:]
  after_breakpoint_ch2 = ch2[breakpoint:]

  new_ch1 = ch1[:breakpoint] + after_breakpoint_ch2
  new_ch2 = ch2[:breakpoint] + after_breakpoint_ch1

  return (new_ch1, new_ch2)

def get_crossover_message(ch1_tuple: Tuple[int, List[int]], ch2_tuple: Tuple[int, List[int]], breakpoint, result1: List[int], result2: List[int]):
  ch1 = "".join(list(map(str, ch1_tuple[1])))
  ch2 = "".join(list(map(str, ch2_tuple[1])))
  result1 = "".join(list(map(str, result1)))
  
  if result2 != None:
    result2 = "".join(list(map(str, result2)))

  ch_len = len(ch1_tuple[1])
  max_padding_right = len(str(ch_len))

  ch1_msg = f'chromosome {ch1_tuple[0] + 1: <{max_padding_right}}: {ch1}\n'
  ch2_msg = f'chromosome {ch2_tuple[0] + 1: <{max_padding_right}}: {ch2}\n'
  delimiter = "-" * ch_len + " breakpoint: {}\n".format(breakpoint)
  res1_msg = (len(f'chromosome {ch1_tuple[0] + 1: <{max_padding_right}}') + 2) * ' ' + result1 + '\n'

  res2_msg = ''
  if result2 != None:
    res2_msg = (len(f'chromosome {ch1_tuple[0] + 1: <{max_padding_right}}') + 2) * ' ' + result2

  return ch1_msg + ch2_msg + delimiter + res1_msg + res2_msg


def perform_crossover(in_crossover: List[Tuple[int, int]], out_crossover: List[Tuple[int, int]], out_file: TextIOWrapper = None) -> List[List[int]]:
  # Taking the first tuple which is of type (index, chromosome).
  chromosome_length = len(in_crossover[0][1])
  
  result = []

  is_special_case = len(in_crossover) == 3
  if is_special_case:
    first, second, third = in_crossover

    breakpoint_1 = random.randint(1, chromosome_length - 1)
    descendant_1 = chromosomes_crossover(first[1], second[1], breakpoint_1)[0]
    result.append((first[0], descendant_1))
    
    breakpoint_2 = random.randint(1, chromosome_length - 1)
    descendant_2 = chromosomes_crossover(second[1], third[1], breakpoint_2)[0]
    result.append((second[0], descendant_2))

    breakpoint_3 = random.randint(1, chromosome_length - 1)
    descendant_3 = chromosomes_crossover(third[1], first[1], breakpoint_3)[0]
    result.append((third[0], descendant_3))

    if out_file != None:
      out_file.write(get_crossover_message(first, second, breakpoint_1, descendant_1, None))
      out_file.write("\n\n")

      out_file.write(get_crossover_message(first, second, breakpoint_2, descendant_2, None))
      out_file.write("\n\n")

      out_file.write(get_crossover_message(first, second, breakpoint_3, descendant_3, None))
      out_file.write("\n\n")

  elif len(in_crossover) % 2 == 1:
    out_crossover.append(in_crossover.pop())

  if is_special_case == False:
    for idx in range(0, len(in_crossover), 2):
      ch1_tuple = in_crossover[idx]
      ch2_tuple = in_crossover[idx + 1]
      breakpoint = random.randint(1, chromosome_length - 1)

      (new_ch1, new_ch2) = chromosomes_crossover(ch1_tuple[1], ch2_tuple[1], breakpoint)
      result.append((ch1_tuple[0], new_ch1))
      result.append((ch2_tuple[0], new_ch2))

      if out_file != None:
        out_file.write(get_crossover_message(ch1_tuple, ch2_tuple, breakpoint, new_ch1, new_ch2))
        out_file.write("\n\n")
  
  result += out_crossover
  
  # Rearranging the items so that they're in order(e.g. from 1 to N).
  temp_res = [] + result
  for (real_idx, ch) in temp_res:
    result[real_idx] = ch

  return result

def perform_mutation(population: List[List[int]], prob_mutation, mutation_type: MutationTypes, out_file: TextIOWrapper = None):
  ch_len = len(population[0])

  if mutation_type == MutationTypes.RARE:
    for idx in range(0, len(population)):
      prob_chromosome = random.uniform(0, 1)
      
      if out_file != None:
        out_file.write("\t")
        out_file.write("chromosome {}: prob_chromosome = {}".format(idx + 1, prob_chromosome))
      
      if prob_chromosome >= prob_mutation:
        if out_file != None:
          out_file.write("\n")
        
        continue
      
      chromosome = population[idx]
      old_chromosome = chromosome[:]
      gene_idx = random.randint(0, ch_len - 1)
      chromosome[gene_idx] = 1 - chromosome[gene_idx]

      if out_file != None:
        out_file.write(" < {} ===> mutation occurred! {} ---idx = {}---> {}\n".format(prob_mutation, "".join(list(map(str, old_chromosome))), gene_idx, "".join(list(map(str, chromosome)))))
  elif mutation_type == MutationTypes.EACH_GENE:
    for idx in range(0, len(population)):
      chromosome = population[idx]
      old_chromosome = chromosome[:]

      has_mutation = False
      # 'current chromosome'

      if out_file != None:
        out_file.write("\t")
        out_file.write("chromosome {}:\n".format(idx + 1))

      for gene_idx in range(0, len(chromosome)):
        prob_gene = random.uniform(0, 1)
        if prob_gene >= prob_mutation:
          continue
        
        has_mutation = True

        old_gene = chromosome[gene_idx]
        chromosome[gene_idx] = 1 - chromosome[gene_idx]
        if out_file != None:
          out_file.write("\t\tgene {} has changed: {} -> {}\n".format(gene_idx + 1, old_gene, chromosome[gene_idx]))

      if out_file != None:
        if has_mutation:
          out_file.write("\t\tfinal chromosome: {} -> {}\n\n".format("".join(list(map(str, old_chromosome))), "".join(list(map(str, chromosome)))))
        else:
          out_file.write("\tnothing changed\n\n")
  
  return population

def get_elitist_chromosome(config, population: List[List[int]]):
  func = config.function

  ch_values = list(map(lambda ch: quadratic_fn(func.a, func.b, func.c, get_chromosome_value(ch, func.lb, func.ub)), population))

  (idx, max_fitness) = reduce(lambda acc, crt: (crt[0] if crt[1] == max(acc[1], crt[1]) else acc[0] , max(acc[1], crt[1])), zip(range(0, len(population)), ch_values))

  return (idx, population[idx], max_fitness)

def get_average_fitness(config, population: List[List[int]]):
  func = config.function

  ch_values = list(map(lambda ch: quadratic_fn(func.a, func.b, func.c, get_chromosome_value(ch, func.lb, func.ub)), population))

  return (reduce(lambda acc, crt: acc + crt, ch_values)) / len(ch_values)

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
