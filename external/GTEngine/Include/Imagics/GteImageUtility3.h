// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Imagics/GteImage3.h>

namespace gte
{

class GTE_IMPEXP ImageUtility3
{
public:
    // All but the Draw* functions are operations on binary images.  Let the
    // image have d0 columns, d1 rows, and d2 slices.  The input image must
    // have zeros on its boundaries x = 0, x = d0-1, y = 0, y = d1-1, z = 0,
    // and z = d2-1.  The 0-valued voxels are considered to be background.
    // The 1-valued voxels are considered to be foreground.  In some of the
    // operations, to save memory and time the input image is modified by the
    // algorithms.  If you need to preserve the input image, make a copy of it
    // before calling these functions.

    // Compute the 6-connected components of a binary image.  The input image
    // is modified to avoid the cost of making a copy.  On output, the image
    // values are the labels for the components.  The array components[k],
    // k >= 1, contains the indices for the k-th component.
    static void GetComponents6(Image3<int>& image,
        std::vector<std::vector<size_t>>& components);

    // Compute the 18-connected components of a binary image.  The input image
    // is modified to avoid the cost of making a copy.  On output, the image
    // values are the labels for the components.  The array components[k],
    // k >= 1, contains the indices for the k-th component.
    static void GetComponents18(Image3<int>& image,
        std::vector<std::vector<size_t>>& components);

    // Compute the 26-connected components of a binary image.  The input image
    // is modified to avoid the cost of making a copy.  On output, the image
    // values are the labels for the components.  The array components[k],
    // k >= 1, contains the indices for the k-th component.
    static void GetComponents26(Image3<int>& image,
        std::vector<std::vector<size_t>>& components);

    // Dilate the image using a structuring element that contains the
    // 6-connected neighbors.
    static void Dilate6(Image3<int> const& inImage, Image3<int>& outImage);

    // Dilate the image using a structuring element that contains the
    // 18-connected neighbors.
    static void Dilate18(Image3<int> const& inImage, Image3<int>& outImage);

    // Dilate the image using a structuring element that contains the
    // 26-connected neighbors.
    static void Dilate26(Image3<int> const& inImage, Image3<int>& outImage);

    // Compute coordinate-directional convex set.  For a given coordinate
    // direction (x, y, or z), identify the first and last 1-valued voxels
    // on a segment of voxels in that direction.  All voxels from first to
    // last are set to 1.  This is done for all segments in each of the
    // coordinate directions.
    static void ComputeCDConvex(Image3<int>& image);

    // Use a depth-first search for filling a 6-connected region.  This is
    // nonrecursive, simulated by using a heap-allocated "stack".  The input
    // (x,y,z) is the seed point that starts the fill.
    template <typename PixelType>
    static void FloodFill6(Image3<PixelType>& image, int x, int y, int z,
        PixelType foreColor, PixelType backColor);

    // Visit pixels using Bresenham's line drawing algorithm.  The callback
    // represents the action you want applied to each voxel as it is visited.
    static void DrawLine(int x0, int y0, int z0, int x1, int y1, int z1,
        std::function<void(int, int, int)> const& callback);

private:
    // Dilation using the specified structuring element.
    static void Dilate(int numNeighbors, std::array<int, 3> const* delta,
        Image3<int> const& inImage, Image3<int>& outImage);

    // Connected component labeling using depth-first search.
    static void GetComponents(int numNeighbors, int const* delta,
        Image3<int>& image, std::vector<std::vector<size_t>>& components);
};


template <typename PixelType>
void ImageUtility3::FloodFill6(Image3<PixelType>& image, int x, int y,
    int z, PixelType foreColor, PixelType backColor)
{
    // Test for a valid seed.
    int const dim0 = image.GetDimension(0);
    int const dim1 = image.GetDimension(1);
    int const dim2 = image.GetDimension(2);
    if (x < 0 || x >= dim0 || y < 0 || y >= dim1 || z < 0 || z >= dim2)
    {
        // The seed point is outside the image domain, so nothing to fill.
        return;
    }

    // Allocate the maximum amount of space needed for the stack.  An empty
    // stack has top == -1.
    size_t const numPixels = image.GetNumPixels();
    int* xStack = new int[numPixels];
    int* yStack = new int[numPixels];
    int* zStack = new int[numPixels];

    // Push seed point onto stack if it has the background color.  All points
    // pushed onto stack have background color backColor.
    int top = 0;
    xStack[top] = x;
    yStack[top] = y;
    zStack[top] = z;

    while (top >= 0)  // stack is not empty
    {
        // Read top of stack.  Do not pop since we need to return to this
        // top value later to restart the fill in a different direction.
        x = xStack[top];
        y = yStack[top];
        z = zStack[top];

        // Fill the pixel.
        image(x, y, z) = foreColor;

        int xp1 = x + 1;
        if (xp1 < dim0 && image(xp1, y, z) == backColor)
        {
            // Push pixel with background color.
            ++top;
            xStack[top] = xp1;
            yStack[top] = y;
            zStack[top] = z;
            continue;
        }

        int xm1 = x - 1;
        if (0 <= xm1 && image(xm1, y, z) == backColor)
        {
            // Push pixel with background color.
            ++top;
            xStack[top] = xm1;
            yStack[top] = y;
            zStack[top] = z;
            continue;
        }

        int yp1 = y + 1;
        if (yp1 < dim1 && image(x, yp1, z) == backColor)
        {
            // Push pixel with background color.
            ++top;
            xStack[top] = x;
            yStack[top] = yp1;
            zStack[top] = z;
            continue;
        }

        int ym1 = y - 1;
        if (0 <= ym1 && image(x, ym1, z) == backColor)
        {
            // Push pixel with background color.
            ++top;
            xStack[top] = x;
            yStack[top] = ym1;
            zStack[top] = z;
            continue;
        }

        int zp1 = z + 1;
        if (zp1 < dim2 && image(x, y, zp1) == backColor)
        {
            // Push pixel with background color.
            ++top;
            xStack[top] = x;
            yStack[top] = y;
            zStack[top] = zp1;
            continue;
        }

        int zm1 = z - 1;
        if (0 <= zm1 && image(x, y, zm1) == backColor)
        {
            // Push pixel with background color.
            ++top;
            xStack[top] = x;
            yStack[top] = y;
            zStack[top] = zm1;
            continue;
        }

        // Done in all directions, pop and return to search other directions.
        --top;
    }

    delete[] xStack;
    delete[] yStack;
    delete[] zStack;
}


}
