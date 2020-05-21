#include "cameras/projective_camera.h"


ProjectiveCamera::ProjectiveCamera(const AnimatedTransform& CameraToWorld,
    const Transform& CameraToScreen, 
    const Eigen::Vector3d& screenWindow, 
    double shutterOpen, 
    double shutterClose, 
    double lensr, 
    double focald, 
    Film* film, 
    const Medium* medium)
{
    // TODO: 374
    /** ScreenToRaster = Scale(film->fullResolution.x,
        film->fullResolution.y, 1) *
        Scale(1 / (screenWindow.pMax.x - screenWindow.pMin.x),
            1 / (screenWindow.pMin.y - screenWindow.pMax.y), 1) *
        Translate(Vector3f(-screenWindow.pMin.x, -screenWindow.pMax.y, 0));
    RasterToScreen = Inverse(ScreenToRaster);
    RasterToCamera = Inverse(CameraToScreen) * RasterToScreen;**/
}
