module SequenceSummary

export normalizeDNA,
       composition,
       gc_content,
       complement,
       reverse_complement,
       parse_fasta,
       kmer_match

# # uncomment the following line if you intend to use BioSequences types
# using BioSequences

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a String.
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    for base in seq
        # note: `N` indicates an unknown base
        occursin(base, "AGCTN") || error("invalid base $base")
    end
    return seq # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end

"""
    composition((::AbstractString)

Creates a dictionary of nucleotide bases (including N)
Uses dictionary to count number of each base in a 
sequence.
Returns dictionary with updated counts for each base
"""

function composition(sequence)
    #sequence = normalizeDNA(sequence) 
    base_comp = Dict()
    #initialize base keys to 0
    for base in "ATCGNS" #S corresponds to G or C base
        base_comp[base] = 0
    end
    #add 1 to value of base key if encountered
    for base in sequence
        if haskey(base_comp, base)
            base_comp[base] += 1
        end    
    end
    
    return base_comp
end

"""
    gc_content(::AbstractString)

Uses the composition function to obtain counts for 
each base (including N) in a sequence.
Returns the GC content by adding counts for G and C and
dividing by the number of bases in the sequence
"""

function gc_content(seq)
    # obtain number of each base in the sequence
    baseComp = composition(seq)
    #calculate GC content
    return (baseComp['G'] + baseComp['C'] + baseComp['S'])/length(seq)
end 

"""
    complement(::AbstractString)

Creates a dictionary where each base is a key with its value 
assigned to the complement base. The complement bases
are compiled to form the sequence complement.
Returns sequence complement
"""
    
function complement(sequence)
    #make a dictionary for complement of each base
        complements = Dict("A" => 'T',
                           "T" => 'A',
                           "G" => 'C',
                           "N" => 'N',
                           "C" => 'G')
    #create the complement of the sequence   
    seq_complement = collect("")
    for base in sequence
        base = uppercase(string(base))
        push!(seq_complement, complements[base])
    end
    seq_complement = join(seq_complement)
    return seq_complement
 end

 """
    reverse_complement(::AbstractString)

Returns the complement of the sequence in reverse 
 """

function reverse_complement(sequence)
    return reverse(complement(sequence))
end

"""
    parse_fasta(path)

Given the path to a file, it separates the contents 
into strings of headers and sequences.
Returns an Array for headers and for sequences
"""
function parse_fasta(path)
    headers = [] #holds tuples of parsed headers
    sequences = [] #holds one long concatendated string/sequence
    seq_holder = [] #temporary holder of current sequence strings
        for line in eachline(path)
            if startswith(line, '>')
                if length(headers) == 0  #if first header is encountered
                    push!(headers, line[2:end]) #add header to headers as a string
                elseif length(headers) > 0 #indicating a new header has been encountered
                    seq_holder = join(seq_holder) #join lines of previous sequence into 1 string
                    push!(sequences, seq_holder)  #add complete string to sequences
                    seq_holder = []               #reset to empty state
                    push!(headers, line[2:end]) #add new header line to headers
                end
            else
                #normalizeDNA(line) #check if only valid bases are present: update- removed temporarily
                push!(seq_holder, line) 
            end
        end
        seq_holder = join(seq_holder)
        push!(sequences, seq_holder)
        return (headers, sequences)  
    end
"""
    kmer_match(::AbstractString, ::Int64)

Using a sequence and provided kmer length, creates a dictionary with
each unique kmer as a key. Returns an array of all of the kmer keys
"""
function kmer_match(seq, k)
    1 <= k <= length(seq) || error("k must be a positive integer less than the length of the sequence")
    kmers = Dict() # initialize dictionary 
    stopindex = (length(seq)) - (k-1) 
    for i in 1:stopindex
        kmer = seq[i:i+(k-1)] # change to index the seq from i to i+k-1
        bases_present = composition(kmer) #find which bases are present in kmer
        #check if kmer consists of only AGCT
        if (bases_present['A'] + bases_present['T'] + bases_present['C'] + bases_present['G'])/length(kmer) == 1
            kmers[kmer] = 1 # create a key for it with a value of 1
        end
    end
    #store and return all of the kmer keys
    all_kmers = [] 
    for (key, value) in kmers
        push!(all_kmers, key)
    end
    return all_kmers
end

"""
    kmer_dist(StringArray, StringArray)

Compares the kmers found in both sequences with the total 
number of kmers in each to determine and return distance metric. 
Distance metric is a value between 0 and 1 where identical 
sequences have a distance metric of 0.
"""
function kmer_dist(set1, set2)
    shared_kmers = intersect(set1, set2)
    kmer_tot = union(set1, set2)
    set_dist = 1 - (length(shared_kmers)/length(kmer_tot))
    return set_dist
end

end # module SequenceSummary