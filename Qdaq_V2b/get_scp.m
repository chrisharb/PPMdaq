function [scope,sn,sn_c,nscp,msg] = get_scp(library)

import LibTiePie.Const.* %Import LibTiePie constants library
import LibTiePie.Enum.* %Import LibTiePie enumeration library
library.DeviceList.update() %Update the Device list
allowedProductIDs = [PID.HS3, PID.HS4, PID.HS4D];
import LibTiePie.Oscilloscope; %Import the LibTiePie Oscilloscope library for creation of a blank instrument
scps = Oscilloscope.empty; %Create an empty Oscilloscope
sn = [];
if library.DeviceList.Count == 0 %Check if MATLAB detected any devices
    scope = -1;
    sn = -1;
    sn_c = -1;
    nscp = -1;
    msg = -1;
elseif ismac %If computer runs Apple IOS this code will not work :-(
    scope = -1;
    sn = -1;
    nscp = -1;
    sn_c = -1;
    msg = -1;
else %Run the code to import the instruments
    for k = 0 : library.DeviceList.Count - 1
        item = library.DeviceList.getItemByIndex(k);
        %if ismember(item.ProductId, allowedProductIDs) && item.canOpen(DEVICETYPE.OSCILLOSCOPE)
            scps(end+1) = item.openOscilloscope();
            sn(end+1) = item.SerialNumber;
            %disp([num2str(scps(end).Interfaces+1) ' channel(s) available on scope '...
            %num2str(k+1) ', model: ' scps(end).Name ', serial number: ' num2str(sn(end))])
        %end
    end
    clear item
    nscp = length(scps);
    if length(scps) > 1
        %disp('Attempting to combine instruments together');
        scope = library.DeviceList.createAndOpenCombinedDevice(scps);
        sn_c = scope.SerialNumber;
        % Remove HS3/HS4(D) objects, not required anymore:
        clear scps;
        msg=['I found ' num2str(length(scope.Channels)) ' available channels across ' num2str(nscp) ' oscilloscopes, which are now combined.'];
    else
        scope = scps;
        sn_c = scope.SerialNumber;
        clear scps;
        msg=['I found ' num2str(length(scope.Channels)) ' available channels which is a standalone instrument'];
    end
end
clear item
