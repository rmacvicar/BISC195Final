module SequenceSummary

export normalizeDNA,
       composition,
       gc_content,
       complement,
       reverse_complement,
       parse_fasta

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
    sequence = normalizeDNA(sequence) 
    base_comp = Dict()
    #initialize base keys to 0
    for base in "ATCGN"
        base_comp[base] = 0
    end
    #add 1 to value of base key if encountered
    for base in sequence
        base_comp[base] += 1     
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
    return (baseComp['G'] + baseComp['C'])/length(seq)
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
                normalizeDNA(line) #check if only valid bases are present
                push!(seq_holder, line) 
            end
        end
        seq_holder = join(seq_holder)
        push!(sequences, seq_holder)
        return (headers, sequences)  
    end
        
end # module SequenceSummary