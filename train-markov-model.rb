#!/usr/bin/env ruby
# Outputs a transition probability matrix for an input FASTA file.
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
               'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
abort("Usage: #{$0} fasta-file output-file") if ARGV.size < 2
# Get a single sequence from a FASTA file, returning nil if EOF
def get_sequence(file)
  sequence = ""
  while line = file.gets
    break if line[0, 1] == '>'
    sequence += line.strip
  end
  sequence = nil if sequence == ""
  sequence
end

infile = File.open(ARGV[0])
infile.gets # remove the 1st defline to make it easier to get sequences, hacky
outfile = File.open(ARGV[1], "w")
transition_counts = {}
while seq = get_sequence(infile)
  unless seq =~ /^[ACDEFGHIKLMNPQRSTVWY]+$/
    warn "Skipping non-standard/wildcard sequence #{seq}"
    next
  end
  prev_char = '0'
  seq.each_char do |c|
    transition_counts[prev_char + '|' + c] ||= 0
    transition_counts[prev_char + '|' + c] += 1
    prev_char = c
  end
end

totals = {}

transition_counts.each { |k, v| totals[k[0,1]] ||= 0; totals[k[0,1]] += v }

transition_probs = transition_counts.clone
transition_probs.each { |k, v| transition_probs[k] = v.to_f/totals[k[0,1]] }

amino_acids.each do |i|
  amino_acids.each do |j|
    prob = transition_probs[i + '|' + j]
    prob = 0.0 unless prob
    outfile << "#{prob}\n"
  end
end
