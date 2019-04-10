#! /usr/bin/env python

"""
Observes the next num_chars characters from a given file handle.
Returns the characters (if possible) and returns the handle back.
"""
def peek(fp, num_chars):
  data = fp.read(num_chars)
  if len(data) == 0:
    return ''
  fp.seek(num_chars * -1, 1)
  return data

"""
Returns a single read from the given FASTA/FASTQ file.
Parameter header contains only the header of the read.
Parameter lines contains all lines of the read, which include:
- header
- seq
- '+' if FASTQ
- quals if FASTQ
Parameter lines is an array of strings, each for one component.
Please note that multiline FASTA/FASTQ entries (e.g. sequence line)
will be truncated into one single line.
Author: Ivan Sovic, 2015.
"""
def get_single_read(fp):
  lines = []

  STATE_HEADER = 0;
  STATE_SEQ = 1;
  STATE_QUAL_SEPARATOR = 2;
  STATE_QUAL = 4;
  state = STATE_HEADER;     # State machine. States:
                            # 0 header, 1 seq, 2 '+' line, 3 quals.
  num_lines = 0;
  header = '';
  seq = '';
  qual_separator = '';
  qual = '';
  lines = [];

  next_char = peek(fp, 1);
  while (len(next_char) > 0):
    line = fp.readline().rstrip();
    next_char = peek(fp, 1);

    if (state == STATE_HEADER):
      if (len(line) == 0): continue;
      header_separator = line[0];
      header = line[1:];			# Strip the '>' or '@' sign from the beginning.
      lines.append(line);
      next_state = STATE_SEQ;
    elif (state == STATE_SEQ):
      seq += line;
      if (len(next_char) == 0):
        lines.append(seq);
        next_state = STATE_HEADER;
        break;      # EOF.
      elif (header_separator == '>' and next_char == header_separator):
        lines.append(seq);
        next_state = STATE_HEADER;
        break;      # This function reads only one sequence.
      elif (header_separator == '@' and next_char == '+'):
        lines.append(seq);
        next_state = STATE_QUAL_SEPARATOR;
      else:
        next_state = STATE_SEQ;

    elif (state == STATE_QUAL_SEPARATOR):
      qual_separator = line;
      lines.append(line);
      next_state = STATE_QUAL;

    elif (state == STATE_QUAL):
      qual += line;
      if (len(next_char) == 0):
        lines.append(qual);
        next_state = STATE_HEADER;
        break;      # EOF.
      elif (next_char == header_separator and len(qual) == len(seq)):
        lines.append(qual);
        next_state = STATE_HEADER;
        break;      # This function reads only one sequence.
      else:
        next_state = STATE_QUAL;
    state = next_state;

  return [header, lines];

def yield_seq(fofn_lines):
  """
  Yields a single sequence from a set of FASTA/FASTQ files provided
  by the list of file names.
  """
  for file_name in fofn_lines:
    fp_in = open(file_name, 'r');

    while(True):
      [header, seq] = get_single_read(fp_in);
      if (len(seq) == 0): break;
      yield(seq)

  fp_in.close();

def main():                 # pragma: no cover
  pass                      # pragma: no cover

if __name__ == "__main__":  # pragma: no cover
  main()                    # pragma: no cover
