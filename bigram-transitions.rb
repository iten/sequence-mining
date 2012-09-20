#!/usr/bin/env ruby
# Outputs a transition probability matrix for input files.
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
               'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
abort("Usage: #{$0} laplace-smooth-amount input-files output-file") if ARGV.size < 3
smooth_amount = ARGV[0].to_f
infilenames = ARGV[1..-2]
outfile = File.open(ARGV.last, "w")
transition_counts = {}
infilenames.each do |filename|
  text = File.open(filename).read
  text.gsub!(/<.*?>/, "")
  text_words = text.split
  text_words.each do |word|
    unless word =~ /^[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]+$/
      next
    end
    word.upcase!
    prev_char = '0'
    word.each_char do |c|
      transition_counts[prev_char + '|' + c] ||= 0
      transition_counts[prev_char + '|' + c] += 1
      prev_char = c
    end
  end
end

totals = {}

# Count # of total transitions for each "prev" amino acid -> all "next" amino
# acids. This is so we can normalize later.
transition_counts.each do |trans, count|
  totals[trans[0,1]] ||= 0.0
  totals[trans[0,1]] += count
end

# Normalize the probability and do Laplace smoothing if desired
transition_probs = transition_counts.clone
transition_counts.each do |trans, count|
  # Laplace smoothing if smooth_amount > 0
  transition_probs[trans] = (count.to_f + smooth_amount) / (totals[trans[0,1]] + 20*smooth_amount)
end

amino_acids.each do |i|
  amino_acids.each do |j|
    prob = transition_probs[i + '|' + j]
    prob = 0.0 unless prob
    outfile << "#{prob}\n"
  end
end
