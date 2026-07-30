function result = verLessThan(package_name, requested_version)
% VERLESSTHAN Minimal Octave compatibility for QETLAB's MATLAB parser guard.
%
% QETLAB asks only whether the MATLAB release predates 8.3 before choosing
% addParameter versus the retired addParamValue spelling. Modern Octave
% supports addParameter but its built-in verLessThan has no "matlab" package.

if ~strcmpi(package_name, 'matlab')
    error('QuantumEntanglementTools:OctaveCompatibility', ...
        'The oracle shim supports only QETLAB''s MATLAB version query.');
end
if ~strcmp(requested_version, '8.3')
    error('QuantumEntanglementTools:OctaveCompatibility', ...
        'Unexpected MATLAB version query: %s.', requested_version);
end
result = false;
end
