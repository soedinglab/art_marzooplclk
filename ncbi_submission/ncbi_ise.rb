#!/usr/bin/env ruby
#
# Clean up a transcriptome sequences to conform with NCBI requirements
#
# Usage: ./ncbi_ise.rb transcriptome.fasta > cleaned_transcriptome.fasta
#

require 'bio'
require 'optparse'
require 'set'

options = {:report => nil, :fcsreport => nil}

parser = OptionParser.new do|opts|
	opts.banner = "Usage: ncbi_ise.rb [options] input.fasta"
	opts.on('-r','--report file','Report file generated when uploading a .sqn to TSA during submission') do |repfile|
		options[:report] = repfile;
	end

	opts.on('-f', '--fcs-report file', 'Report file generated after a TSA submission is processed but still has errors') do |repfile|
		options[:fcsreport] = repfile;
	end

	opts.on('-h', '--help', 'Displays Help') do
		puts opts
		exit
	end
end

parser.parse!

if ARGV[0].nil?
	$stderr.write "You must supply an input file\n" 
    puts parser
    exit
end

input_fasta = Bio::FlatFile.auto(ARGV[0])

# Find out the lengths of each sequence
#
sequence_lengths = {}
input_fasta.each do |entry|  
	eid=entry.entry_id
	sequence_lengths[eid] = entry.naseq.length
end

input_fasta = Bio::FlatFile.auto(ARGV[0])


# Parse an ncbi error report file to find sequences for trimming
#
#
# For example we are parsing lines like this
#
# ERROR fhd.sqn FhD07614 35..35 VECTOR_MATCH File: fhd.sqn, Code(VECTOR_MATCH), Sequence-id: FhD07614, Interval: 35..35, ...
#
# And our aim is to extra the id FhD07614 and trim interval 35..35
#
trim_sequences = {}
unless ( options[:report].nil? )
	File.open(options[:report], "r").each do |line|
		if line !~ /Sequence-id: (.*), Interval: (.*),/
			$stderr.write "Unrecognised VecScreen output \n#{line}"
		else
			id,interval = line.match(/Sequence-id: (.*), Interval: (.*),/).captures
			vstart,vend = interval.split("..").collect { |e| e.to_i-1 }		
			trim_sequences[id] = [vstart,vend]
		end
	end
end

# Parse an ncbi FSC report to find more sequences for trimming or deletion
#
# For example we are parsing line like this for trimming
#
# FhD01419	1316	75..203	vector/etc
#
#
duplicated_sequences = Set.new
excluded_sequences = Set.new
unless ( options[:fcsreport].nil? )
	File.open(options[:fcsreport], "r").each do |line|  
		parts = line.split("\t")
		
		# Match and parse lines specifying parts of sequences to trim
		#
		if (parts.length == 4) && ( parts[2] =~ /\d+\.\.\d+/)
			id,interval = parts[0],parts[2]
			vstart,vend = interval.split("..").collect { |e| e.to_i-1 }		
			trim_sequences[id] = [vstart,vend]
		end
		
		# Match lines specifying duplicated sequences
		if (line =~ /^lcl\|/ )
			# Parse the ids of duplicates and them to the duplicates list
			dups = line.to_enum(:scan, /lcl\|([^~]*)/).map { Regexp.last_match.captures[0] }
			dlens = dups.collect { |dp| sequence_lengths[dp] }
			rep = dups[dlens.index(dlens.max)]
			dups.delete(rep)
			duplicated_sequences.merge(dups)
		end
		
		# This is designed to match lines like this
		# comp34504_c0_seq1	206	vector/etc
		#
		if ( parts.length == 3) && ( parts[1] =~ /^\d+$/ )
			excluded_sequences.add(parts[0])	
		end

	end
end


too_short=0
toomany_n=0
trim_n=0
dup_n = 0
excluded_n = 0

input_fasta.each do |entry|  

	eid=entry.entry_id


	# No X's or *'s
	#
	seq = entry.naseq.gsub(/[Xx]/,"n")
	seq.sub!(/^n*/,"")
	seq.sub!(/n*$/,"")

	skip=false

	if trim_sequences.has_key?(eid)
		seq.slice!(trim_sequences[eid][0],trim_sequences[eid][1])
		trim_n+=1		
	end

	# require 'byebug'; byebug

	# Skip sequences that are too short or have too many N's
	if (seq.length < 200)
		too_short+=1
		skip=true
	elsif seq =~ /nnnnnnnnnnnnnn/ || 
		seq =~ /NNNNNNNNNNNNNN/ ||    # No more than 14 N's anywhere
		seq[1..20] =~ /NNNNNNNNNN/ || 
		seq[1..20] =~ /nnnnnnnnnn/ || # No more than 10 N's at start or end
		seq[-20..-1] =~ /NNNNNNNNNN/ || 
		seq[-20..-1] =~ /nnnnnnnnnn/ || 
		(seq.count('n')/seq.length) > 0.1 || 
		(seq.count('N')/seq.length) > 0.1
		toomany_n+=1
		skip=true
	end

	# Skip duplicates
	if duplicated_sequences.member?(eid)
		dup_n += 1
		skip=true
	end
	
	# Skip excludes
	if excluded_sequences.member?(eid)
		excluded_n += 1
		skip=true
	end

	$stdout.write ">#{entry.entry_id}\n#{seq.upcase}\n" unless skip



end

$stderr.write "Skipped #{too_short} sequences because they were too short\n"
$stderr.write "Skipped #{toomany_n} sequences because they had too many Ns\n"
$stderr.write "Trimmed #{trim_n} sequences\n"
$stderr.write "Skipped #{dup_n} sequences with too much redundancy\n"
$stderr.write "Skipped #{excluded_n} sequences identified for exclusion\n"
