function meta = parse_lab_identifiers(path_and_file)
% Parses Run/Lot/Sample/Waveguide/Date identifiers out of a lab file path,
% using the folder-token and filename regex conventions used across the
% "on Samples" test data tree (e.g. .../Run G/AM156/SD51/id4.1/...).
%
%   meta = parse_lab_identifiers(path_and_file)
%
% meta fields (all string): Run, Lot, Sample, Waveguide, Date
% Unresolved fields fall back to "UnknownRun" / "UnknownLot" / "UnknownSample"
% / "UnknownWG"; Date falls back to "" if no date-like token is found.
%
% Promoted, unmodified, from parseAllIdentifiers() in
% OverheadLoss_IntermediateData_DefectDetection.m so both that script and
% compile_overhead_images.m can share one robust implementation instead of
% fragile hardcoded path-token indexing.

    meta.Run       = "UnknownRun";
    meta.Lot       = "UnknownLot";
    meta.Sample    = "UnknownSample";
    meta.Waveguide = "UnknownWG";
    meta.Date      = "";

    parts = split(path_and_file, filesep);
    fname = parts{end};

    for k = 1:numel(parts)
        token = parts{k};

        if regexp(token,'^(Run|GL)\w*','once')
            meta.Run = token;
        end
        if regexp(token,'^AM\d+','once')
            meta.Lot = token;
        end
        if regexp(token,'^SD\d+','once')
            meta.Sample = token;
        end
        if regexp(token,'^(WG\d+|id\d+(\.\d+)?)','once')
            meta.Waveguide = token;
        end
    end

    if meta.Waveguide == "UnknownWG"
        m = regexp(fname,'(WG\d+|id\d+(\.\d+)?)','match','once');
        if ~isempty(m), meta.Waveguide = m; end
    end

    d = regexp(path_and_file,'\d{4}[-_]\d{2}[-_]\d{2}','match','once');
    if ~isempty(d)
        meta.Date = strrep(d,'_','-');
    end

    meta.Waveguide = normalizeWaveguide(meta.Waveguide);

end

function wg = normalizeWaveguide(wg)
    wg = string(wg);
    if startsWith(wg,"WG")
        n = regexp(wg,'\d+','match','once');
        wg = "WG" + string(str2double(n));
    end
end
