#!/usr/bin/env ruby
# Output to stdout a visualization of a transition matrix
puts "Usage: #{$0} matrix [precision]" if ARGV.size == 0
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
               'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
s = File.open(ARGV[0]).read
if ARGV.size >= 2
  precision = ARGV[1].to_i
else
  precision = 4
end
probs = s.split("\n")
# Print top row of labels
amino_acids.each do |aa|
  print " " * (precision+2) + "#{aa}"
end
print "\n"
# Print each row, with label
amino_acids.each_with_index do |aa_prev, i|
  print "#{aa_prev} "
  amino_acids.each_with_index do |aa_next, j|
    print " %.#{precision}f" % probs[i*20 + j].to_f
  end
  print "\n"
end
