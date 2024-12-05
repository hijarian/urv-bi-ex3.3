% FINGERPRINT MATCHING SCORE
% Argument:   ff1 -  First transformed Fingerprint 
%                  data format: array of minutiae, each one is a row [x y theta class]
%             ff2 -  Second transformed Fingerprint
%                  data format: array of minutiae, each one is a row [x y theta class]
% Returns:    Distance
function Distance = Match_fingerprints(ff1, ff2)
    MaxD = 0;
    TotalDistance = 0;

    ff1Size = size(ff1, 1);
    ff2Size = size(ff2, 1);

    isMatchedFF1 = false(size(ff1, 1), 1); % Logical array to mark elements as matched
    isMatchedFF2 = false(size(ff2, 1), 1); % Logical array to mark elements as matched

    % Loop through each minutia in ff1
    for i = 1:ff1Size
        if all(isMatchedFF2)
            % no more minutiae in ff2 to match with
            break;
        end

        mm1 = ff1(i, :);

        % temporary variables to store the best match between mm1 and mm2
        % to pick if no "exact" match was found by `mm`.
        CurrentMinDistance = inf;
        BestMatchIndex = -1;

        % Loop through each minutia in ff2
        for j = 1:ff2Size
            if isMatchedFF2(j)
                % skip already matched minutiae
                continue;
            end

            mm2 = ff2(j, :);
            % actually perform the matching
            [IsCloseEnough, CandidateDistance] = mm(mm1, mm2);

            if CandidateDistance > MaxD
                % Max distance is a global value across all the possible combinations so immediately update it
                MaxD = CandidateDistance;
            end

            if IsCloseEnough
                if CandidateDistance < CurrentMinDistance
                    CurrentMinDistance = CandidateDistance;
                    BestMatchIndex = j;
                end
            end % if IsMatch
        end % for j

        if BestMatchIndex ~= -1
            TotalDistance = TotalDistance + CurrentMinDistance;
            isMatchedFF2(BestMatchIndex) = true; % Mark as matched
            isMatchedFF1(i) = true; % Mark as matched
        end
    end

    % Calculate the number of leftover minutiae
    nm1 = sum(~isMatchedFF1);
    nm2 = sum(~isMatchedFF2);

    % Calculate the final distance
    K = max(nm1, nm2) * MaxD;
    TotalSize = ff1Size + ff2Size;
    Distance = TotalDistance + K / TotalSize;
end

% each minutia is a row in the matrix
% each row has 4 columns: x, y, theta, class
function SpatialDistance = sd(LeftMinutia, RightMinutia)
    dx = LeftMinutia(1) - RightMinutia(1);
    dy = LeftMinutia(2) - RightMinutia(2);
    SpatialDistance = sqrt(dx^2 + dy^2);
end

% each minutia is a row in the matrix
% each row has 4 columns: x, y, theta, class
% NOTE: IMPORTANT: this is ULTRA dependent on the possible values of theta!!!
% THIS MUST BE IN SYNC WITH THE FUNCTION GETTING THE THETA VALUES
function DirectionalDistance = dd(LeftMinutia, RightMinutia)
    DTheta = abs(LeftMinutia(4) - RightMinutia(4));
    DTheta = min(DTheta, 2*pi - DTheta);
    DirectionalDistance = DTheta;
end

% each minutia is a row in the matrix
% each row has 4 columns: x, y, theta, class
% note that CLASS IS THE LAST COLUMN
function [IsMatch, Distance] = mm(LeftMinutia, RightMinutia)
    DISTANCE_TOLERANCE = 10; % hand-crafted based on the input data ranges
    THETA_TOLERANCE = 1; % hand-crafted based on the input data ranges
    SD_WEIGHT = 1; % arbitrarily chosen
    DD_WEIGHT = 1; % arbitrarily chosen

    SpatialDistance = sd(LeftMinutia, RightMinutia);
    DirectionalDistance = dd(LeftMinutia, RightMinutia);
    IsSameType = LeftMinutia(4) == RightMinutia(4);

    IsMatch = SpatialDistance < DISTANCE_TOLERANCE && DirectionalDistance < THETA_TOLERANCE && IsSameType;
    Distance = SD_WEIGHT * SpatialDistance + DD_WEIGHT * DirectionalDistance;
end
