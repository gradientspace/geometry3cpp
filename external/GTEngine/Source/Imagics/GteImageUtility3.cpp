// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2018
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <Imagics/GteImageUtility3.h>
using namespace gte;

void ImageUtility3::GetComponents6(Image3<int>& image,
    std::vector<std::vector<size_t>>& components)
{
    std::array<int, 6> neighbors;
    image.GetNeighborhood(neighbors);
    GetComponents(6, &neighbors[0], image, components);
}

void ImageUtility3::GetComponents18(Image3<int>& image,
    std::vector<std::vector<size_t>>& components)
{
    std::array<int, 18> neighbors;
    image.GetNeighborhood(neighbors);
    GetComponents(18, &neighbors[0], image, components);
}

void ImageUtility3::GetComponents26(Image3<int>& image,
    std::vector<std::vector<size_t>>& components)
{
    std::array<int, 26> neighbors;
    image.GetNeighborhood(neighbors);
    GetComponents(26, &neighbors[0], image, components);
}

void ImageUtility3::Dilate6(Image3<int> const& inImage, Image3<int>& outImage)
{
    std::array<std::array<int, 3>, 6> neighbors;
    inImage.GetNeighborhood(neighbors);
    Dilate(6, &neighbors[0], inImage, outImage);
}

void ImageUtility3::Dilate18(Image3<int> const& inImage,
    Image3<int>& outImage)
{
    std::array<std::array<int, 3>, 18> neighbors;
    inImage.GetNeighborhood(neighbors);
    Dilate(18, &neighbors[0], inImage, outImage);
}

void ImageUtility3::Dilate26(Image3<int> const& inImage,
    Image3<int>& outImage)
{
    std::array<std::array<int, 3>, 26> neighbors;
    inImage.GetNeighborhood(neighbors);
    Dilate(26, &neighbors[0], inImage, outImage);
}

void ImageUtility3::ComputeCDConvex(Image3<int>& image)
{
    int const dim0 = image.GetDimension(0);
    int const dim1 = image.GetDimension(1);
    int const dim2 = image.GetDimension(2);

    Image3<int> temp = image;
    int i0, i1, i2;
    for (i1 = 0; i1 < dim1; ++i1)
    {
        for (i0 = 0; i0 < dim0; ++i0)
        {
            int i2min;
            for (i2min = 0; i2min < dim2; ++i2min)
            {
                if ((temp(i0, i1, i2min) & 1) == 0)
                {
                    temp(i0, i1, i2min) |= 2;
                }
                else
                {
                    break;
                }
            }
            if (i2min < dim2)
            {
                int i2max;
                for (i2max = dim2 - 1; i2max >= i2min; --i2max)
                {
                    if ((temp(i0, i1, i2max) & 1) == 0)
                    {
                        temp(i0, i1, i2max) |= 2;
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }
    }

    for (i2 = 0; i2 < dim2; ++i2)
    {
        for (i0 = 0; i0 < dim0; ++i0)
        {
            int i1min;
            for (i1min = 0; i1min < dim1; ++i1min)
            {
                if ((temp(i0, i1min, i2) & 1) == 0)
                {
                    temp(i0, i1min, i2) |= 2;
                }
                else
                {
                    break;
                }
            }
            if (i1min < dim1)
            {
                int i1max;
                for (i1max = dim1 - 1; i1max >= i1min; --i1max)
                {
                    if ((temp(i0, i1max, i2) & 1) == 0)
                    {
                        temp(i0, i1max, i2) |= 2;
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }
    }

    for (i2 = 0; i2 < dim2; ++i2)
    {
        for (i1 = 0; i1 < dim1; ++i1)
        {
            int i0min;
            for (i0min = 0; i0min < dim0; ++i0min)
            {
                if ((temp(i0min, i1, i2) & 1) == 0)
                {
                    temp(i0min, i1, i2) |= 2;
                }
                else
                {
                    break;
                }
            }
            if (i0min < dim0)
            {
                int i0max;
                for (i0max = dim0 - 1; i0max >= i0min; --i0max)
                {
                    if ((temp(i0max, i1, i2) & 1) == 0)
                    {
                        temp(i0max, i1, i2) |= 2;
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }
    }

    for (size_t i = 0; i < image.GetNumPixels(); ++i)
    {
        image[i] = (temp[i] & 2 ? 0 : 1);
    }
}

void ImageUtility3::DrawLine(int x0, int y0, int z0, int x1, int y1, int z1,
    std::function<void(int, int, int)> const& callback)
{
    // Starting point of line.
    int x = x0, y = y0, z = z0;

    // Direction of line.
    int dx = x1-x0, dy = y1-y0, dz = z1-z0;

    // Increment or decrement depending on direction of line.
    int sx = (dx > 0 ? 1 : (dx < 0 ? -1 : 0));
    int sy = (dy > 0 ? 1 : (dy < 0 ? -1 : 0));
    int sz = (dz > 0 ? 1 : (dz < 0 ? -1 : 0));

    // Decision parameters for voxel selection.
    if (dx < 0)
    {
        dx = -dx;
    }
    if (dy < 0)
    {
        dy = -dy;
    }
    if (dz < 0)
    {
        dz = -dz;
    }
    int ax = 2*dx, ay = 2*dy, az = 2*dz;
    int decX, decY, decZ;

    // Determine largest direction component, single-step related variable.
    int maxValue = dx, var = 0;
    if (dy > maxValue)
    {
        maxValue = dy;
        var = 1;
    }
    if (dz > maxValue)
    {
        var = 2;
    }

    // Traverse Bresenham line.
    switch (var)
    {
    case 0:  // Single-step in x-direction.
        decY = ay - dx;
        decZ = az - dx;
        for (/**/; /**/; x += sx, decY += ay, decZ += az)
        {
            // Process voxel.
            callback(x, y, z);

            // Take Bresenham step.
            if (x == x1)
            {
                break;
            }
            if (decY >= 0)
            {
                decY -= ax;
                y += sy;
            }
            if (decZ >= 0)
            {
                decZ -= ax;
                z += sz;
            }
        }
        break;
    case 1:  // Single-step in y-direction.
        decX = ax - dy;
        decZ = az - dy;
        for (/**/; /**/; y += sy, decX += ax, decZ += az)
        {
            // Process voxel.
            callback(x, y, z);

            // Take Bresenham step.
            if (y == y1)
            {
                break;
            }
            if (decX >= 0)
            {
                decX -= ay;
                x += sx;
            }
            if (decZ >= 0)
            {
                decZ -= ay;
                z += sz;
            }
        }
        break;
    case 2:  // Single-step in z-direction.
        decX = ax - dz;
        decY = ay - dz;
        for (/**/; /**/; z += sz, decX += ax, decY += ay)
        {
            // Process voxel.
            callback(x, y, z);

            // Take Bresenham step.
            if (z == z1)
            {
                break;
            }
            if (decX >= 0)
            {
                decX -= az;
                x += sx;
            }
            if (decY >= 0)
            {
                decY -= az;
                y += sy;
            }
        }
        break;
    }
}

void ImageUtility3::Dilate(int numNeighbors, std::array<int, 3> const* delta,
    Image3<int> const& inImage, Image3<int>& outImage)
{
    int const bound0M1 = inImage.GetDimension(0) - 1;
    int const bound1M1 = inImage.GetDimension(1) - 1;
    int const bound2M1 = inImage.GetDimension(2) - 1;
    for (int i2 = 1; i2 < bound2M1; ++i2)
    {
        for (int i1 = 1; i1 < bound1M1; ++i1)
        {
            for (int i0 = 1; i0 < bound0M1; ++i0)
            {
                if (inImage(i0, i1, i2) == 0)
                {
                    for (int n = 0; n < numNeighbors; ++n)
                    {
                        int d0 = delta[n][0];
                        int d1 = delta[n][1];
                        int d2 = delta[n][2];
                        if (inImage(i0 + d0, i1 + d1, i2 + d2) == 1)
                        {
                            outImage(i0, i1, i2) = 1;
                            break;
                        }
                    }
                }
                else
                {
                    outImage(i0, i1, i2) = 1;
                }
            }
        }
    }
}

void ImageUtility3::GetComponents(int numNeighbors, int const* delta,
    Image3<int>& image, std::vector<std::vector<size_t>>& components)
{
    size_t const numVoxels = image.GetNumPixels();
    std::vector<int> numElements(numVoxels);
    std::vector<size_t> vstack(numVoxels);
    size_t i, numComponents = 0;
    int label = 2;
    for (i = 0; i < numVoxels; ++i)
    {
        if (image[i] == 1)
        {
            int top = -1;
            vstack[++top] = i;

            int& count = numElements[numComponents + 1];
            count = 0;
            while (top >= 0)
            {
                size_t v = vstack[top];
                image[v] = -1;
                int j;
                for (j = 0; j < numNeighbors; ++j)
                {
                    size_t adj = v + delta[j];
                    if (image[adj] == 1)
                    {
                        vstack[++top] = adj;
                        break;
                    }
                }
                if (j == numNeighbors)
                {
                    image[v] = label;
                    ++count;
                    --top;
                }
            }

            ++numComponents;
            ++label;
        }
    }

    if (numComponents > 0)
    {
        components.resize(numComponents + 1);
        for (i = 1; i <= numComponents; ++i)
        {
            components[i].resize(numElements[i]);
            numElements[i] = 0;
        }

        for (i = 0; i < numVoxels; ++i)
        {
            int value = image[i];
            if (value != 0)
            {
                // Labels started at 2 to support the depth-first search,
                // so they need to be decremented for the correct labels.
                image[i] = --value;
                components[value][numElements[value]] = i;
                ++numElements[value];
            }
        }
    }
}
