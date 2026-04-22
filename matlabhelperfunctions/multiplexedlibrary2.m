function seqs = multiplexedlibrary2(n, l, m)
% Generate DNA index sequences using Illumina-style design principles.
%   n : number of sequences
%   l : length of each sequence (>= 4)
%   m : minimum edit distance (Levenshtein) between any two sequences
%
% The function aims to:
%   - keep per-position base balance (A/T/C/G ~ 1:1:1:1, as close as possible)
%   - avoid long homopolymer runs (e.g. AAA, TTT)
%   - enforce a minimum pairwise edit distance m

    % --- 1. Input validation ---
    if l < 4
        error('Sequence length l must be at least 4.');
    end
    if m >= l
        warning('Minimum edit distance m is >= sequence length l. This may be very hard or impossible to satisfy.');
    end

    alphabet = 'ATCG';                % order: A, T, C, G
    numAlphabet = numel(alphabet);
    
    % store results as an n x l char array
    seqs = repmat('N', n, l);  
    
    % cumulative base counts per position: 4 x l (rows correspond to A/T/C/G)
    base_counts = zeros(numAlphabet, l); 
    
    maxTrials = 1e6;                  % safety limit to avoid infinite loops

    fprintf('--- Index generation started (n=%d, l=%d, m=%d) ---\n', n, l, m);

    for i = 1:n
        trial = 0;
        while true
            trial = trial + 1;
            if trial > maxTrials
                error(['Could not find a sequence satisfying all constraints. ', ...
                       'Try reducing m or increasing l.']);
            end
            
            cand = repmat('N', 1, l);
            
            % --- 2. Generate one candidate sequence of length l ---
            for j = 1:l
                % Base counts at this position so far
                current_counts = base_counts(:, j);
                
                % Start with all bases as candidate pool
                available_bases_idx = 1:numAlphabet;
                
                % B. Homopolymer constraint: avoid triples like AAA, TTT, etc.
                if j >= 3 
                    last_base        = cand(j-1);
                    second_last_base = cand(j-2);
                    
                    if last_base == second_last_base
                        % If the previous two bases are identical (e.g. 'CC'),
                        % disallow that base at the current position.
                        homo_base_idx = find(alphabet == last_base);
                        available_bases_idx = setdiff(available_bases_idx, homo_base_idx);
                        
                        if isempty(available_bases_idx)
                            % If homopolymer constraint leaves no base,
                            % fall back to any base except the homopolymer one.
                            temp_pool = setdiff(1:numAlphabet, homo_base_idx); 
                            selected_idx = temp_pool(randi(numel(temp_pool)));
                            cand(j) = alphabet(selected_idx);
                            continue; % move on to next position j
                        end
                    end
                end
                
                % C. Balance constraint (greedy):
                %    Among the allowed bases, choose one with minimal usage
                %    count at this position so far to keep A/T/C/G balanced.
                
                if isempty(available_bases_idx)
                    % Safety fallback (should rarely happen)
                    selected_idx = randi(numAlphabet);
                    cand(j) = alphabet(selected_idx);
                else
                    current_counts_in_pool = current_counts(available_bases_idx);
                    min_count = min(current_counts_in_pool);
                    
                    % all bases in the pool that share the minimal count
                    min_bases_rel_idx = find(current_counts_in_pool == min_count);
                    
                    % map back to indices in the full alphabet
                    selected_idx = available_bases_idx( ...
                        min_bases_rel_idx(randi(numel(min_bases_rel_idx))) );
                    cand(j) = alphabet(selected_idx);
                end
            end % for j=1:l
            
            % --- 3. Check EDIT distance (Levenshtein) to all previous sequences ---
            is_dist_ok = true;
            if i > 1
                prevSeqs = seqs(1:i-1, :);
                for k = 1:i-1
                    d = editDistance(prevSeqs(k, :), cand);
                    if d < m
                        is_dist_ok = false;
                        break;
                    end
                end
            end
            
            % --- 4. Accept candidate and update counts if it passes ---
            if is_dist_ok
                seqs(i, :) = cand;
                
                % Update base_counts for this accepted sequence
                for j = 1:l
                   base_index = find(alphabet == cand(j));
                   base_counts(base_index, j) = base_counts(base_index, j) + 1;
                end
                
                break;
            end
        end % while (trial)
    end % for (i=1:n)
    
    % --- Final report ---
    fprintf('\n--- Summary of generated library ---\n');
    disp(['Requested minimum edit distance (m): ', num2str(m)]);
    disp('Final base counts per position (rows = A/T/C/G):');
    disp(base_counts);
    
    % quantify how balanced each position is: max - min count at that column
    per_pos_spread = max(base_counts, [], 1) - min(base_counts, [], 1);
    fprintf('Per-position base count spread (max - min):\n');
    disp(per_pos_spread);
    fprintf('Maximum spread over all positions: %d\n', max(per_pos_spread));
    fprintf('Smaller spread means closer to perfect A/T/C/G balance at each position.\n');
end

% -------------------------------------------------------------------------
% Local helper: Levenshtein edit distance between two strings
% -------------------------------------------------------------------------
function d = editDistance(s1, s2)
    s1 = char(s1);
    s2 = char(s2);
    n1 = length(s1);
    n2 = length(s2);
    
    % DP matrix (n1+1) x (n2+1)
    dp = zeros(n1+1, n2+1);
    
    % initialization
    dp(:,1) = 0:n1;  % cost of deletions
    dp(1,:) = 0:n2;  % cost of insertions
    
    for i = 2:n1+1
        for j = 2:n2+1
            cost = (s1(i-1) ~= s2(j-1)); % 0 if same char, 1 if different
            dp(i,j) = min([ ...
                dp(i-1,j)   + 1,  ... % deletion
                dp(i,j-1)   + 1,  ... % insertion
                dp(i-1,j-1) + cost ...% substitution
            ]);
        end
    end
    
    d = dp(n1+1, n2+1);
end