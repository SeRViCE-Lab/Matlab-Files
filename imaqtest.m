flag.Visualize = true;

%cleanupWorkspace

obj.colorVideoDevice = imaq.VideoDevice('kinect', 1);
set(obj.colorVideoDevice, 'ReturnDataTyoe', 'uint8');
obj.depthVideoDevice = imaq.VideoDevice('kinect', 2);