#define NM_TO_A 10
#define CAL2JOULE 4.184

extern "C" __device__ void imageParticles(const real4 * box, const real3 * referenceParticle, real3 * particleToImage)
{
    // Periodic boundary conditions imaging of particleToImage with respect to referenceParticle

    real3 delta = *referenceParticle - *particleToImage;

    real3 offset;
    offset.x = floor(delta.x / (*box).x + 0.5f) * (*box).x;
    offset.y = floor(delta.y / (*box).y + 0.5f) * (*box).y;
    offset.z = floor(delta.z / (*box).z + 0.5f) * (*box).z;

    *particleToImage += offset;
}
