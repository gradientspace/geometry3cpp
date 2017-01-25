// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <Imagics/GteImage2.h>

namespace gte
{

class GTE_IMPEXP ImageUtility2
{
public:
    // All but the Draw* functions are operations on binary images.  Let the
    // image have d0 columns and d1 rows.  The input image must have zeros on
    // its boundaries x = 0, x = d0-1, y = 0, and y = d1-1.  The 0-valued
    // pixels are considered to be background.  The 1-valued pixels are
    // considered to be foreground.  In some of the operations, to save memory
    // and time the input image is modified by the algorithms.  If you need
    // to preserve the input image, make a copy of it before calling these
    // functions.
    //
    // Dilation and erosion functions do not have the requirement that the
    // boundary pixels of the binary image inputs be zero.

    // Compute the 4-connected components of a binary image.  The input image
    // is modified to avoid the cost of making a copy.  On output, the image
    // values are the labels for the components.  The array components[k],
    // k >= 1, contains the indices for the k-th component.
    static void GetComponents4(Image2<int>& image,
        std::vector<std::vector<size_t>>& components);

    // Compute the 8-connected components of a binary image.  The input image
    // is modified to avoid the cost of making a copy.  On output, the image
    // values are the labels for the components.  The array components[k],
    // k >= 1, contains the indices for the k-th component.
    static void GetComponents8(Image2<int>& image,
        std::vector<std::vector<size_t>>& components);

    // Compute a dilation with a structuring element consisting of the
    // 4-connected neighbors of each pixel.  The input image is binary with 0
    // for background and 1 for foreground.  The output image must be an
    // object different from the input image.
    static void Dilate4(Image2<int> const& input, Image2<int>& output);

    // Compute a dilation with a structuring element consisting of the
    // 8-connected neighbors of each pixel.  The input image is binary with 0
    // for background and 1 for foreground.  The output image must be an
    // object different from the input image.
    static void Dilate8(Image2<int> const& input, Image2<int>& output);

    // Compute a dilation with a structing element consisting of neighbors
    // specified by offsets relative to the pixel.  The input image is binary
    // with 0 for background and 1 for foreground.  The output image must be
    // an object different from the input image.
    static void Dilate(Image2<int> const& input, int numNeighbors,
        std::array<int, 2> const* neighbors, Image2<int>& output);

    // Compute an erosion with a structuring element consisting of the
    // 4-connected neighbors of each pixel.  The input image is binary with 0
    // for background and 1 for foreground.  The output image must be an
    // object different from the input image.  If zeroExterior is true, the
    // image exterior is assumed to be 0, so 1-valued boundary pixels are set
    // to 0; otherwise, boundary pixels are set to 0 only when they have
    // neighboring image pixels that are 0.
    static void Erode4(Image2<int> const& input, bool zeroExterior,
        Image2<int>& output);

    // Compute an erosion with a structuring element consisting of the
    // 8-connected neighbors of each pixel.  The input image is binary with 0
    // for background and 1 for foreground.  The output image must be an
    // object different from the input image.  If zeroExterior is true, the
    // image exterior is assumed to be 0, so 1-valued boundary pixels are
    // set to 0; otherwise, boundary pixels are set to 0 only when they have
    // neighboring image pixels that are 0.
    static void Erode8(Image2<int> const& input, bool zeroExterior,
        Image2<int>& output);

    // Compute an erosion with a structuring element consisting of neighbors
    // specified by offsets relative to the pixel.  The input image is binary
    // with 0 for background and 1 for foreground.  The output image must be
    // an object different from the input image.  If zeroExterior is true, the
    // image exterior is assumed to be 0, so 1-valued boundary pixels are set
    // to 0; otherwise, boundary pixels are set to 0 only when they have
    // neighboring image pixels that are 0.
    static void Erode(Image2<int> const& input, bool zeroExterior,
        int numNeighbors, std::array<int, 2> const* neighbors,
        Image2<int>& output);

    // Compute an opening with a structuring element consisting of the
    // 4-connected neighbors of each pixel.  The input image is binary with 0
    // for background and 1 for foreground.  The output image must be an
    // object different from the input image.  If zeroExterior is true, the
    // image exterior is assumed to consist of 0-valued pixels; otherwise,
    // the image exterior is assumed to consist of 1-valued pixels.
    static void Open4(Image2<int> const& input, bool zeroExterior,
        Image2<int>& output);

    // Compute an opening with a structuring element consisting of the
    // 8-connected neighbors of each pixel.  The input image is binary with 0
    // for background and 1 for foreground.  The output image must be an
    // object different from the input image.  If zeroExterior is true, the
    // image exterior is assumed to consist of 0-valued pixels; otherwise,
    // the image exterior is assumed to consist of 1-valued pixels.
    static void Open8(Image2<int> const& input, bool zeroExterior,
        Image2<int>& output);

    // Compute an opening with a structuring element consisting of neighbors
    // specified by offsets relative to the pixel.  The input image is binary
    // with 0 for background and 1 for foreground.  The output image must be
    // an object different from the input image.  If zeroExterior is true, the
    // image exterior is assumed to consist of 0-valued pixels; otherwise,
    // the image exterior is assumed to consist of 1-valued pixels.
    static void Open(Image2<int> const& input, bool zeroExterior,
        int numNeighbors, std::array<int, 2> const* neighbors,
        Image2<int>& output);

    // Compute a closing with a structuring element consisting of the
    // 4-connected neighbors of each pixel.  The input image is binary with 0
    // for background and 1 for foreground.  The output image must be an
    // object different from the input image.  If zeroExterior is true, the
    // image exterior is assumed to consist of 0-valued pixels; otherwise,
    // the image exterior is assumed to consist of 1-valued pixels.
    static void Close4(Image2<int> const& input, bool zeroExterior,
        Image2<int>& output);

    // Compute a closing with a structuring element consisting of the
    // 8-connected neighbors of each pixel.  The input image is binary with 0
    // for background and 1 for foreground.  The output image must be an
    // object different from the input image.  If zeroExterior is true, the
    // image exterior is assumed to consist of 0-valued pixels; otherwise,
    // the image exterior is assumed to consist of 1-valued pixels.
    static void Close8(Image2<int> const& input, bool zeroExterior,
        Image2<int>& output);

    // Compute a closing with a structuring element consisting of neighbors
    // specified by offsets relative to the pixel.  The input image is binary
    // with 0 for background and 1 for foreground.  The output image must be
    // an object different from the input image.  If zeroExterior is true, the
    // image exterior is assumed to consist of 0-valued pixels; otherwise,
    // the image exterior is assumed to consist of 1-valued pixels.
    static void Close(Image2<int> const& input, bool zeroExterior,
        int numNeighbors, std::array<int, 2> const* neighbors,
        Image2<int>& output);

    // Locate a pixel and walk around the edge of a component.  The input
    // (x,y) is where the search starts for a nonzero pixel.  If (x,y) is
    // outside the component, the walk is around the outside the component.
    // If the component has a hole and (x,y) is inside that hole, the walk
    // is around the boundary surrounding the hole.  The function returns
    // 'true' on a success walk.  The return value is 'false' when no
    // boundary was found from the starting (x,y).
    static bool ExtractBoundary(int x, int y, Image2<int>& image,
        std::vector<size_t>& boundary);

    // Use a depth-first search for filling a 4-connected region.  This is
    // nonrecursive, simulated by using a heap-allocated "stack".  The input
    // (x,y) is the seed point that starts the fill.
    template <typename PixelType>
    static void FloodFill4(Image2<PixelType>& image, int x, int y,
        PixelType foreColor, PixelType backColor);

    // Compute the L1-distance transform of the binary image.  The function
    // returns the maximum distance and a point at which the maximum
    // distance is attained.
    static void GetL1Distance(Image2<int>& image, int& maxDistance,
        int& xMax, int& yMax);

    // Compute the L2-distance transform of the binary image.  The maximum
    // distance should not be larger than 100, so you have to ensure this is
    // the case for the input image.  The function returns the maximum
    // distance and a point at which the maximum distance is attained.
    // Comments about the algorithm are in the source file.
    static void GetL2Distance(Image2<int> const& image, float& maxDistance,
        int& xMax, int& yMax, Image2<float>& transform);

    // Compute a skeleton of a binary image.  Boundary pixels are trimmed from
    // the object one layer at a time based on their adjacency to interior
    // pixels.  At each step the connectivity and cycles of the object are
    // preserved.  The skeleton overwrites the contents of the input image.
    static void GetSkeleton(Image2<int>& image);

    // In the remaining public member functions, the callback represents the
    // action you want applied to each pixel as it is visited.

    // Visit pixels in a (2*thick+1)x(2*thick+1) square centered at (x,y).
    static void DrawThickPixel(int x, int y, int thick,
        std::function<void(int, int)> const& callback);

    // Visit pixels using Bresenham's line drawing algorithm.
    static void DrawLine(int x0, int y0, int x1, int y1,
        std::function<void(int, int)> const& callback);

    // Visit pixels using Bresenham's circle drawing algorithm.  Set 'solid'
    // to false for drawing only the circle.  Set 'solid' to true to draw all
    // pixels on and inside the circle.
    static void DrawCircle(int xCenter, int yCenter, int radius, bool solid,
        std::function<void(int, int)> const& callback);

    // Visit pixels in a rectangle of the specified dimensions.  Set 'solid'
    // to false for drawing only the rectangle.  Set 'solid' to true to draw
    // all pixels on and inside the rectangle.
    static void DrawRectangle(int xMin, int yMin, int xMax, int yMax,
        bool solid, std::function<void(int, int)> const& callback);

    // Use a depth-first search for filling a 4-connected region.  This is
    // nonrecursive, simulated by using a heap-allocated "stack".  The input
    // (x,y) is the seed point that starts the fill.  The x-value is in
    // {0..xSize-1} and the y-value is in {0..ySize-1}.
    template <typename PixelType>
    static void DrawFloodFill4(int x, int y, int xSize, int ySize,
        PixelType foreColor, PixelType backColor,
        std::function<void(int, int, PixelType)> const& setCallback,
        std::function<PixelType(int, int)> const& getCallback);

private:
    // Connected component labeling using depth-first search.
    static void GetComponents(int numNeighbors, int const* delta,
        Image2<int>& image, std::vector<std::vector<size_t>>& components);

    // Support for GetL2Distance.
    static void L2Check(int x, int y, int dx, int dy, Image2<int>& xNear,
        Image2<int>& yNear, Image2<int>& dist);

    // Support for GetSkeleton.
    static bool Interior2 (Image2<int>& image, int x, int y);
    static bool Interior3 (Image2<int>& image, int x, int y);
    static bool Interior4 (Image2<int>& image, int x, int y);
    static bool MarkInterior (Image2<int>& image, int value,
        bool (*function)(Image2<int>&,int,int));
    static bool IsArticulation (Image2<int>& image, int x, int y);
    static bool ClearInteriorAdjacent (Image2<int>& image, int value);
    static int const msArticulation[256];
};


template <typename PixelType>
void ImageUtility2::FloodFill4(Image2<PixelType>& image, int x, int y,
    PixelType foreColor, PixelType backColor)
{
    // Test for a valid seed.
    int const dim0 = image.GetDimension(0);
    int const dim1 = image.GetDimension(1);
    if (x < 0 || x >= dim0 || y < 0 || y >= dim1)
    {
        // The seed point is outside the image domain, so nothing to fill.
        return;
    }

    // Allocate the maximum amount of space needed for the stack.  An empty
    // stack has top == -1.
    size_t const numPixels = image.GetNumPixels();
    int* xStack = new int[numPixels];
    int* yStack = new int[numPixels];

    // Push seed point onto stack if it has the background color.  All points
    // pushed onto stack have background color backColor.
    int top = 0;
    xStack[top] = x;
    yStack[top] = y;

    while (top >= 0)  // stack is not empty
    {
        // Read top of stack.  Do not pop since we need to return to this
        // top value later to restart the fill in a different direction.
        x = xStack[top];
        y = yStack[top];

        // Fill the pixel.
        image(x, y) = foreColor;

        int xp1 = x + 1;
        if (xp1 < dim0 && image(xp1, y) == backColor)
        {
            // Push pixel with background color.
            ++top;
            xStack[top] = xp1;
            yStack[top] = y;
            continue;
        }

        int xm1 = x - 1;
        if (0 <= xm1 && image(xm1, y) == backColor)
        {
            // Push pixel with background color.
            ++top;
            xStack[top] = xm1;
            yStack[top] = y;
            continue;
        }

        int yp1 = y + 1;
        if (yp1 < dim1 && image(x, yp1) == backColor)
        {
            // Push pixel with background color.
            ++top;
            xStack[top] = x;
            yStack[top] = yp1;
            continue;
        }

        int ym1 = y - 1;
        if (0 <= ym1 && image(x, ym1) == backColor)
        {
            // Push pixel with background color.
            ++top;
            xStack[top] = x;
            yStack[top] = ym1;
            continue;
        }

        // Done in all directions, pop and return to search other directions.
        --top;
    }

    delete[] xStack;
    delete[] yStack;
}

template <typename PixelType>
void ImageUtility2::DrawFloodFill4(int x, int y, int xSize, int ySize,
    PixelType foreColor, PixelType backColor,
    std::function<void(int, int, PixelType)> const& setCallback,
    std::function<PixelType(int, int)> const& getCallback)
{
    // Test for a valid seed.
    if (x < 0 || x >= xSize || y < 0 || y >= ySize)
    {
        // The seed point is outside the image domain, so nothing to fill.
        return;
    }

    // Allocate the maximum amount of space needed for the stack.  An empty
    // stack has top == -1.
    int const numPixels = xSize * ySize;
    std::vector<int> xStack(numPixels), yStack(numPixels);

    // Push seed point onto stack if it has the background color.  All points
    // pushed onto stack have background color backColor.
    int top = 0;
    xStack[top] = x;
    yStack[top] = y;

    while (top >= 0)  // stack is not empty
    {
        // Read top of stack.  Do not pop since we need to return to this
        // top value later to restart the fill in a different direction.
        x = xStack[top];
        y = yStack[top];

        // Fill the pixel.
        setCallback(x, y, foreColor);

        int xp1 = x + 1;
        if (xp1 < xSize && getCallback(xp1, y) == backColor)
        {
            // Push pixel with background color.
            ++top;
            xStack[top] = xp1;
            yStack[top] = y;
            continue;
        }

        int xm1 = x - 1;
        if (0 <= xm1 && getCallback(xm1, y) == backColor)
        {
            // Push pixel with background color.
            ++top;
            xStack[top] = xm1;
            yStack[top] = y;
            continue;
        }

        int yp1 = y + 1;
        if (yp1 < ySize && getCallback(x, yp1) == backColor)
        {
            // Push pixel with background color.
            ++top;
            xStack[top] = x;
            yStack[top] = yp1;
            continue;
        }

        int ym1 = y - 1;
        if (0 <= ym1 && getCallback(x, ym1) == backColor)
        {
            // Push pixel with background color.
            ++top;
            xStack[top] = x;
            yStack[top] = ym1;
            continue;
        }

        // Done in all directions, pop and return to search other directions.
        --top;
    }
}


}
