// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#include <GTEnginePCH.h>
#include <Imagics/GteImageUtility2.h>
#include <limits>
using namespace gte;


void ImageUtility2::GetComponents4(Image2<int>& image,
    std::vector<std::vector<size_t>>& components)
{
    std::array<int, 4> neighbors;
    image.GetNeighborhood(neighbors);
    GetComponents(4, &neighbors[0], image, components);
}

void ImageUtility2::GetComponents8(Image2<int>& image,
    std::vector<std::vector<size_t>>& components)
{
    std::array<int, 8> neighbors;
    image.GetNeighborhood(neighbors);
    GetComponents(8, &neighbors[0], image, components);
}

void ImageUtility2::Dilate4(Image2<int> const& input, Image2<int>& output)
{
    std::array<std::array<int, 2>, 4> neighbors;
    input.GetNeighborhood(neighbors);
    Dilate(input, 4, &neighbors[0], output);
}

void ImageUtility2::Dilate8(Image2<int> const& input, Image2<int>& output)
{
    std::array<std::array<int, 2>, 8> neighbors;
    input.GetNeighborhood(neighbors);
    Dilate(input, 8, &neighbors[0], output);
}

void ImageUtility2::Dilate(Image2<int> const& input, int numNeighbors,
    std::array<int, 2> const* neighbors, Image2<int>& output)
{
    // If the assertion is triggered, the function will run but the output
    // will not be correct.
    LogAssert(&output != &input, "Input and output must be different.");

    output = input;

    // If the pixel at (x,y) is 1, then the pixels at (x+dx,y+dy) are set to 1
    // where (dx,dy) is in the 'neighbors' array.  Boundary testing is used to
    // avoid accessing out-of-range pixels.
    int const dim0 = input.GetDimension(0);
    int const dim1 = input.GetDimension(1);
    for (int y = 0; y < dim1; ++y)
    {
        for (int x = 0; x < dim0; ++x)
        {
            if (input(x, y) == 1)
            {
                for (int j = 0; j < numNeighbors; ++j)
                {
                    int xNbr = x + neighbors[j][0];
                    int yNbr = y + neighbors[j][1];
                    if (0 <= xNbr && xNbr < dim0 && 0 <= yNbr && yNbr < dim1)
                    {
                        output(xNbr, yNbr) = 1;
                    }
                }
            }
        }
    }
}

void ImageUtility2::Erode4(Image2<int> const& input, bool zeroExterior,
    Image2<int>& output)
{
    std::array<std::array<int, 2>, 4> neighbors;
    input.GetNeighborhood(neighbors);
    Erode(input, zeroExterior, 4, &neighbors[0], output);
}

void ImageUtility2::Erode8(Image2<int> const& input, bool zeroExterior,
    Image2<int>& output)
{
    std::array<std::array<int, 2>, 8> neighbors;
    input.GetNeighborhood(neighbors);
    Erode(input, zeroExterior, 8, &neighbors[0], output);
}

void ImageUtility2::Erode(Image2<int> const& input, bool zeroExterior,
    int numNeighbors, std::array<int, 2> const* neighbors,
    Image2<int>& output)
{
    // If the assertion is triggered, the function will run but the output
    // will not be correct.
    LogAssert(&output != &input, "Input and output must be different.");

    output = input;

    // If the pixel at (x,y) is 1, it is changed to 0 when at least one
    // neighbor (x+dx,y+dy) is 0, where (dx,dy) is in the 'neighbors'
    // array.
    int const dim0 = input.GetDimension(0);
    int const dim1 = input.GetDimension(1);
    for (int y = 0; y < dim1; ++y)
    {
        for (int x = 0; x < dim0; ++x)
        {
            if (input(x, y) == 1)
            {
                for (int j = 0; j < numNeighbors; ++j)
                {
                    int xNbr = x + neighbors[j][0];
                    int yNbr = y + neighbors[j][1];
                    if (0 <= xNbr && xNbr < dim0 && 0 <= yNbr && yNbr < dim1)
                    {
                        if (input(xNbr, yNbr) == 0)
                        {
                            output(x, y) = 0;
                            break;
                        }
                    }
                    else if (zeroExterior)
                    {
                        output(x, y) = 0;
                        break;
                    }
                }
            }
        }
    }
}

void ImageUtility2::Open4(Image2<int> const& input, bool zeroExterior,
    Image2<int>& output)
{
    Image2<int> temp(input.GetDimension(0), input.GetDimension(1));
    Erode4(input, zeroExterior, temp);
    Dilate4(temp, output);
}

void ImageUtility2::Open8(Image2<int> const& input, bool zeroExterior,
    Image2<int>& output)
{
    Image2<int> temp(input.GetDimension(0), input.GetDimension(1));
    Erode8(input, zeroExterior, temp);
    Dilate8(temp, output);
}

void ImageUtility2::Open(Image2<int> const& input, bool zeroExterior,
    int numNeighbors, std::array<int, 2> const* neighbors,
    Image2<int>& output)
{
    Image2<int> temp(input.GetDimension(0), input.GetDimension(1));
    Erode(input, zeroExterior, numNeighbors, neighbors, temp);
    Dilate(temp, numNeighbors, neighbors, output);
}

void ImageUtility2::Close4(Image2<int> const& input, bool zeroExterior,
    Image2<int>& output)
{
    Image2<int> temp(input.GetDimension(0), input.GetDimension(1));
    Dilate4(input, temp);
    Erode4(temp, zeroExterior, output);
}

void ImageUtility2::Close8(Image2<int> const& input, bool zeroExterior,
    Image2<int>& output)
{
    Image2<int> temp(input.GetDimension(0), input.GetDimension(1));
    Dilate8(input, temp);
    Erode8(temp, zeroExterior, output);
}

void ImageUtility2::Close(Image2<int> const& input, bool zeroExterior,
    int numNeighbors, std::array<int, 2> const* neighbors,
    Image2<int>& output)
{
    Image2<int> temp(input.GetDimension(0), input.GetDimension(1));
    Dilate(input, numNeighbors, neighbors, temp);
    Erode(temp, zeroExterior, numNeighbors, neighbors, output);
}

bool ImageUtility2::ExtractBoundary(int x, int y, Image2<int>& image,
    std::vector<size_t>& boundary)
{
    // Find a first boundary pixel.
    size_t const numPixels = image.GetNumPixels();
    size_t i;
    for (i = image.GetIndex(x, y); i < numPixels; ++i)
    {
        if (image[i])
        {
            break;
        }
    }
    if (i == numPixels)
    {
        // No boundary pixel found.
        return false;
    }

    int const dx[8] = { -1,  0, +1, +1, +1,  0, -1, -1 };
    int const dy[8] = { -1, -1, -1,  0, +1, +1, +1,  0 };

    // Create a new point list that contains the first boundary point.
    boundary.push_back(i);

    // The direction from background 0 to boundary pixel 1 is (dx[7],dy[7]).
    std::array<int,2> coord = image.GetCoordinates(i);
    int x0 = coord[0], y0 = coord[1];
    int cx = x0, cy = y0;
    int nx = x0-1, ny = y0, dir = 7;

    // Traverse the boundary in clockwise order.  Mark visited pixels as 2.
    image(cx, cy) = 2;
    bool notDone = true;
    while (notDone)
    {
        int j, nbr;
        for (j = 0, nbr = dir; j < 8; ++j, nbr = (nbr + 1) % 8)
        {
            nx = cx + dx[nbr];
            ny = cy + dy[nbr];
            if (image(nx, ny))  // next boundary pixel found
            {
                break;
            }
        }

        if (j == 8)  // (cx,cy) is isolated
        {
            notDone = false;
            continue;
        }

        if (nx == x0 && ny == y0)  // boundary traversal completed
        {
            notDone = false;
            continue;
        }

        // (nx,ny) is next boundary point, add point to list.  Note that
        // the index for the pixel is computed for the original image, not
        // for the larger temporary image.
        boundary.push_back(image.GetIndex(nx, ny));

        // Mark visited pixels as 2.
        image(nx, ny) = 2;

        // Start search for next point.
        cx = nx;
        cy = ny;
        dir = (j + 5 + dir) % 8;
    }

    return true;
}

void ImageUtility2::GetL1Distance(Image2<int>& image, int& maxDistance,
    int& xMax, int& yMax)
{
    int const dim0 = image.GetDimension(0);
    int const dim1 = image.GetDimension(1);
    int const dim0m1 = dim0 - 1;
    int const dim1m1 = dim1 - 1;

    // Use a grass-fire approach, computing distance from boundary to
    // interior one pass at a time.
    bool changeMade = true;
    int distance;
    for (distance = 1, xMax = 0, yMax = 0; changeMade; ++distance)
    {
        changeMade = false;
        int distanceP1 = distance + 1;
        for (int y = 1; y < dim1m1; ++y)
        {
            for (int x = 1; x < dim0m1; ++x)
            {
                if (image(x, y) == distance)
                {
                    if (image(x-1, y) >= distance
                    &&  image(x+1, y) >= distance
                    &&  image(x, y-1) >= distance
                    &&  image(x, y+1) >= distance)
                    {
                        image(x, y) = distanceP1;
                        xMax = x;
                        yMax = y;
                        changeMade = true;
                    }
                }
            }
        }
    }

    maxDistance = --distance;
}

void ImageUtility2::GetL2Distance(Image2<int> const& image,
    float& maxDistance, int& xMax, int& yMax, Image2<float>& transform)
{
    // This program calculates the Euclidean distance transform of a binary
    // input image.  The adaptive algorithm is guaranteed to give exact
    // distances for all distances < 100.  Algorithm sent to me by John Gauch.
    //
    // From John Gauch at University of Kansas:
    // The basic idea is similar to a EDT described recently in PAMI by
    // Laymarie from McGill.  By keeping the dx and dy offset to the nearest
    // edge (feature) point in the image, we can search to see which dx dy is
    // closest to a given point by examining a set of neighbors.  The Laymarie
    // method (and Borgfors) look at a fixed 3x3 or 5x5 neighborhood and call
    // it a day.  What we did was calculate (painfully) what neighborhoods you
    // need to look at to guarentee that the exact distance is obtained.  Thus,
    // you will see in the code, that we L2Check the current distance and
    // depending on what we have so far, we extend the search region.  Since
    // our algorithm for L2Checking the exactness of each neighborhood is on
    // the order N^4, we have only gone to N=100.  In theory, you could make
    // this large enough to get all distances exact.  We have implemented the
    // algorithm to get all distances < 100 to be exact. 
    int const dim0 = image.GetDimension(0);
    int const dim1 = image.GetDimension(1);
    int const dim0m1 = dim0 - 1;
    int const dim1m1 = dim1 - 1;
    int x, y, distance;

    // Create and initialize intermediate images.
    Image2<int> xNear(dim0, dim1);
    Image2<int> yNear(dim0, dim1);
    Image2<int> dist(dim0, dim1);
    for (y = 0; y < dim1; ++y)
    {
        for (x = 0; x < dim0; ++x)
        {
            if (image(x, y) != 0)
            {
                xNear(x, y) = 0;
                yNear(x, y) = 0;
                dist(x, y) = std::numeric_limits<int>::max();
            }
            else
            {
                xNear(x, y) = x;
                yNear(x, y) = y;
                dist(x, y) = 0;
            }
        }
    }

    int const K1 = 1;
    int const K2 = 169;   // 13^2
    int const K3 = 961;   // 31^2
    int const K4 = 2401;  // 49^2
    int const K5 = 5184;  // 72^2

    // Pass in the ++ direction.
    for (y = 0; y < dim1; ++y)
    {
        for (x = 0; x < dim0; ++x)
        {
            distance = dist(x, y);
            if (distance > K1)
            { 
                L2Check(x, y, -1,  0, xNear, yNear, dist); 
                L2Check(x, y, -1, -1, xNear, yNear, dist); 
                L2Check(x, y,  0, -1, xNear, yNear, dist); 
            }
            if (distance > K2)
            { 
                L2Check(x, y, -2, -1, xNear, yNear, dist); 
                L2Check(x, y, -1, -2, xNear, yNear, dist); 
            }
            if (distance > K3)
            { 
                L2Check(x, y, -3, -1, xNear, yNear, dist); 
                L2Check(x, y, -3, -2, xNear, yNear, dist); 
                L2Check(x, y, -2, -3, xNear, yNear, dist); 
                L2Check(x, y, -1, -3, xNear, yNear, dist); 
            }
            if (distance > K4)
            { 
                L2Check(x, y, -4, -1, xNear, yNear, dist); 
                L2Check(x, y, -4, -3, xNear, yNear, dist); 
                L2Check(x, y, -3, -4, xNear, yNear, dist); 
                L2Check(x, y, -1, -4, xNear, yNear, dist); 
            }
            if (distance > K5)
            { 
                L2Check(x, y, -5, -1, xNear, yNear, dist); 
                L2Check(x, y, -5, -2, xNear, yNear, dist); 
                L2Check(x, y, -5, -3, xNear, yNear, dist); 
                L2Check(x, y, -5, -4, xNear, yNear, dist);
                L2Check(x, y, -4, -5, xNear, yNear, dist); 
                L2Check(x, y, -2, -5, xNear, yNear, dist); 
                L2Check(x, y, -3, -5, xNear, yNear, dist); 
                L2Check(x, y, -1, -5, xNear, yNear, dist); 
            }
        }
    }

    // Pass in -- direction.
    for (y = dim1m1; y >= 0; --y)
    {
        for (x = dim0m1; x >= 0; --x)
        {
            distance = dist(x, y);
            if (distance > K1)
            { 
                L2Check(x, y, 1, 0, xNear, yNear, dist); 
                L2Check(x, y, 1, 1, xNear, yNear, dist); 
                L2Check(x, y, 0, 1, xNear, yNear, dist); 
            }
            if (distance > K2)
            { 
                L2Check(x, y, 2, 1, xNear, yNear, dist); 
                L2Check(x, y, 1, 2, xNear, yNear, dist); 
            }
            if (distance > K3)
            { 
                L2Check(x, y, 3, 1, xNear, yNear, dist); 
                L2Check(x, y, 3, 2, xNear, yNear, dist); 
                L2Check(x, y, 2, 3, xNear, yNear, dist); 
                L2Check(x, y, 1, 3, xNear, yNear, dist); 
            }
            if (distance > K4)
            { 
                L2Check(x, y, 4, 1, xNear, yNear, dist); 
                L2Check(x, y, 4, 3, xNear, yNear, dist); 
                L2Check(x, y, 3, 4, xNear, yNear, dist); 
                L2Check(x, y, 1, 4, xNear, yNear, dist); 
            }
            if (distance > K5)
            { 
                L2Check(x, y, 5, 1, xNear, yNear, dist); 
                L2Check(x, y, 5, 2, xNear, yNear, dist); 
                L2Check(x, y, 5, 3, xNear, yNear, dist); 
                L2Check(x, y, 5, 4, xNear, yNear, dist);
                L2Check(x, y, 4, 5, xNear, yNear, dist); 
                L2Check(x, y, 2, 5, xNear, yNear, dist); 
                L2Check(x, y, 3, 5, xNear, yNear, dist); 
                L2Check(x, y, 1, 5, xNear, yNear, dist); 
            }
        }
    }

    // Pass in the +- direction.
    for (y = dim1m1; y >= 0; --y)
    {
        for (x = 0; x < dim0; ++x)
        {
            distance = dist(x, y);
            if (distance > K1)
            { 
                L2Check(x, y, -1, 0, xNear, yNear, dist); 
                L2Check(x, y, -1, 1, xNear, yNear, dist); 
                L2Check(x, y,  0, 1, xNear, yNear, dist); 
            }
            if (distance > K2)
            { 
                L2Check(x, y, -2, 1, xNear, yNear, dist); 
                L2Check(x, y, -1, 2, xNear, yNear, dist); 
            }
            if (distance > K3)
            { 
                L2Check(x, y, -3, 1, xNear, yNear, dist); 
                L2Check(x, y, -3, 2, xNear, yNear, dist); 
                L2Check(x, y, -2, 3, xNear, yNear, dist); 
                L2Check(x, y, -1, 3, xNear, yNear, dist); 
            }
            if (distance > K4)
            { 
                L2Check(x, y, -4, 1, xNear, yNear, dist); 
                L2Check(x, y, -4, 3, xNear, yNear, dist); 
                L2Check(x, y, -3, 4, xNear, yNear, dist); 
                L2Check(x, y, -1, 4, xNear, yNear, dist); 
            }
            if (distance > K5)
            { 
                L2Check(x, y, -5, 1, xNear, yNear, dist); 
                L2Check(x, y, -5, 2, xNear, yNear, dist); 
                L2Check(x, y, -5, 3, xNear, yNear, dist); 
                L2Check(x, y, -5, 4, xNear, yNear, dist);
                L2Check(x, y, -4, 5, xNear, yNear, dist); 
                L2Check(x, y, -2, 5, xNear, yNear, dist); 
                L2Check(x, y, -3, 5, xNear, yNear, dist); 
                L2Check(x, y, -1, 5, xNear, yNear, dist); 
            }
        }
    }

    // Pass in the -+ direction.
    for (y = 0; y < dim1; ++y)
    {
        for (x = dim0m1; x >= 0; --x)
        {
            distance = dist(x, y);
            if (distance > K1)
            { 
                L2Check(x, y, 1,  0, xNear, yNear, dist); 
                L2Check(x, y, 1, -1, xNear, yNear, dist); 
                L2Check(x, y, 0, -1, xNear, yNear, dist); 
            }
            if (distance > K2)
            { 
                L2Check(x, y, 2, -1, xNear, yNear, dist); 
                L2Check(x, y, 1, -2, xNear, yNear, dist); 
            }
            if (distance > K3)
            { 
                L2Check(x, y, 3, -1, xNear, yNear, dist); 
                L2Check(x, y, 3, -2, xNear, yNear, dist); 
                L2Check(x, y, 2, -3, xNear, yNear, dist); 
                L2Check(x, y, 1, -3, xNear, yNear, dist); 
            }
            if (distance > K4)
            { 
                L2Check(x, y, 4, -1, xNear, yNear, dist); 
                L2Check(x, y, 4, -3, xNear, yNear, dist); 
                L2Check(x, y, 3, -4, xNear, yNear, dist); 
                L2Check(x, y, 1, -4, xNear, yNear, dist); 
            }
            if (distance > K5)
            { 
                L2Check(x, y, 5, -1, xNear, yNear, dist); 
                L2Check(x, y, 5, -2, xNear, yNear, dist); 
                L2Check(x, y, 5, -3, xNear, yNear, dist); 
                L2Check(x, y, 5, -4, xNear, yNear, dist);
                L2Check(x, y, 4, -5, xNear, yNear, dist); 
                L2Check(x, y, 2, -5, xNear, yNear, dist); 
                L2Check(x, y, 3, -5, xNear, yNear, dist); 
                L2Check(x, y, 1, -5, xNear, yNear, dist); 
            }
        }
    }

    xMax = 0;
    yMax = 0;
    maxDistance = 0.0f;
    for (y = 0; y < dim1; ++y)
    {
        for (x = 0; x < dim0; ++x)
        {
            float fdistance = sqrt((float)dist(x, y));
            if (fdistance > maxDistance)
            {
                maxDistance = fdistance;
                xMax = x;
                yMax = y;
            }
            transform(x, y) = fdistance;
        }
    }
}

void ImageUtility2::GetSkeleton(Image2<int>& image)
{
    int const dim0 = image.GetDimension(0);
    int const dim1 = image.GetDimension(1);

    // Trim pixels, mark interior as 4.
    bool notDone = true;
    while (notDone)
    {
        if (MarkInterior(image, 4, Interior4))
        {
            // No interior pixels, trimmed set is at most 2-pixels thick.
            notDone = false;
            continue;
        }

        if (ClearInteriorAdjacent(image, 4))
        {
            // All remaining interior pixels are either articulation points
            // or part of blobs whose boundary pixels are all articulation
            // points.  An example of the latter case is shown below.  The
            // background pixels are marked with '.' rather than '0' for
            // readability.  The interior pixels are marked with '4' and the
            // boundary pixels are marked with '1'.
            //
            //   .........
            //   .....1...
            //   ..1.1.1..
            //   .1.141...
            //   ..14441..
            //   ..1441.1.
            //   .1.11.1..
            //   ..1..1...
            //   .........
            //
            // This is a pathological problem where there are many small holes
            // (0-pixel with north, south, west, and east neighbors all
            // 1-pixels) that your application can try to avoid by an initial
            // pass over the image to fill in such holes.  Of course, you do
            // have problems with checkerboard patterns...
            notDone = false;
            continue;
        }
    }

    // Trim pixels, mark interior as 3.
    notDone = true;
    while (notDone)
    {
        if (MarkInterior(image, 3, Interior3))
        {
            // No interior pixels, trimmed set is at most 2-pixels thick.
            notDone = false;
            continue;
        }

        if (ClearInteriorAdjacent(image, 3))
        {
            // All remaining 3-values can be safely removed since they are
            // not articulation points and the removal will not cause new
            // holes.
            for (int y = 0; y < dim1; ++y)
            {
                for (int x = 0; x < dim0; ++x)
                {
                    if (image(x, y) == 3 && !IsArticulation(image, x, y))
                    {
                        image(x, y) = 0;
                    }
                }
            }
            notDone = false;
            continue;
        }
    }

    // Trim pixels, mark interior as 2.
    notDone = true;
    while (notDone)
    {
        if (MarkInterior(image, 2, Interior2))
        {
            // No interior pixels, trimmed set is at most 1-pixel thick.
            // Call it a skeleton.
            notDone = false;
            continue;
        }

        if (ClearInteriorAdjacent(image, 2))
        {
            // Removes 2-values that are not articulation points.
            for (int y = 0; y < dim1; ++y)
            {
                for (int x = 0; x < dim0; ++x)
                {
                    if (image(x, y) == 2 && !IsArticulation(image, x, y))
                    {
                        image(x, y) = 0;
                    }
                }
            }
            notDone = false;
            continue;
        }
    }

    // Make the skeleton a binary image.
    size_t const numPixels = image.GetNumPixels();
    for (size_t i = 0; i < numPixels; ++i)
    {
        if (image[i] != 0)
        {
            image[i] = 1;
        }
    }
}

void ImageUtility2::DrawThickPixel(int x, int y, int thick,
    std::function<void(int, int)> const& callback)
{
    for (int dy = -thick; dy <= thick; ++dy)
    {
        for (int dx = -thick; dx <= thick; ++dx)
        {
            callback(x + dx, y + dy);
        }
    }
}

void ImageUtility2::DrawLine(int x0, int y0, int x1, int y1,
    std::function<void(int, int)> const& callback)
{
    // Starting point of line.
    int x = x0, y = y0;

    // Direction of line.
    int dx = x1 - x0, dy = y1 - y0;

    // Increment or decrement depending on direction of line.
    int sx = (dx > 0 ? 1 : (dx < 0 ? -1 : 0));
    int sy = (dy > 0 ? 1 : (dy < 0 ? -1 : 0));

    // Decision parameters for pixel selection.
    if (dx < 0)
    {
        dx = -dx;
    }
    if (dy < 0)
    {
        dy = -dy;
    }
    int ax = 2 * dx, ay = 2 * dy;
    int decX, decY;

    // Determine largest direction component, single-step related variable.
    int maxValue = dx, var = 0;
    if (dy > maxValue)
    {
        var = 1;
    }

    // Traverse Bresenham line.
    switch (var)
    {
    case 0:  // Single-step in x-direction.
        decY = ay - dx;
        for (/**/; /**/; x += sx, decY += ay)
        {
            callback(x, y);

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
        }
        break;
    case 1:  // Single-step in y-direction.
        decX = ax - dy;
        for (/**/; /**/; y += sy, decX += ax)
        {
            callback(x, y);

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
        }
        break;
    }
}

void ImageUtility2::DrawCircle(int xCenter, int yCenter, int radius,
    bool solid, std::function<void(int, int)> const& callback)
{
    int x, y, dec;

    if (solid)
    {
        int xValue, yMin, yMax, i;
        for (x = 0, y = radius, dec = 3 - 2*radius; x <= y; ++x)
        {
            xValue = xCenter + x;
            yMin = yCenter - y;
            yMax = yCenter + y;
            for (i = yMin; i <= yMax; ++i)
            {
                callback(xValue, i);
            }

            xValue = xCenter - x;
            for (i = yMin; i <= yMax; ++i)
            {
                callback(xValue, i);
            }

            xValue = xCenter + y;
            yMin = yCenter - x;
            yMax = yCenter + x;
            for (i = yMin; i <= yMax; ++i)
            {
                callback(xValue, i);
            }

            xValue = xCenter - y;
            for (i = yMin; i <= yMax; ++i)
            {
                callback(xValue, i);
            }

            if (dec >= 0)
            {
                dec += -4*(y--) + 4;
            }
            dec += 4*x + 6;
        }
    }
    else
    {
        for (x = 0, y = radius, dec = 3 - 2*radius; x <= y; ++x)
        {
            callback(xCenter + x, yCenter + y);
            callback(xCenter + x, yCenter - y);
            callback(xCenter - x, yCenter + y);
            callback(xCenter - x, yCenter - y);
            callback(xCenter + y, yCenter + x);
            callback(xCenter + y, yCenter - x);
            callback(xCenter - y, yCenter + x);
            callback(xCenter - y, yCenter - x);

            if (dec >= 0)
            {
                dec += -4*(y--) + 4;
            }
            dec += 4*x + 6;
        }
    }
}

void ImageUtility2::DrawRectangle(int xMin, int yMin, int xMax,
    int yMax, bool solid, std::function<void(int, int)> const& callback)
{
    int x, y;

    if (solid)
    {
        for (y = yMin; y <= yMax; ++y)
        {
            for (x = xMin; x <= xMax; ++x)
            {
                callback(x, y);
            }
        }
    }
    else
    {
        for (x = xMin; x <= xMax; ++x)
        {
            callback(x, yMin);
            callback(x, yMax);
        }
        for (y = yMin + 1; y <= yMax - 1; ++y)
        {
            callback(xMin, y);
            callback(xMax, y);
        }
    }
}

void ImageUtility2::GetComponents(int numNeighbors, int const* delta,
    Image2<int>& image, std::vector<std::vector<size_t>>& components)
{
    size_t const numPixels = image.GetNumPixels();
    int* numElements = new int[numPixels];
    size_t* vstack = new size_t[numPixels];
    size_t i, numComponents = 0;
    int label = 2;
    for (i = 0; i < numPixels; ++i)
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
    delete[] vstack;

    if (numComponents > 0)
    {
        components.resize(numComponents + 1);
        for (i = 1; i <= numComponents; ++i)
        {
            components[i].resize(numElements[i]);
            numElements[i] = 0;
        }

        for (i = 0; i < numPixels; ++i)
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
    delete[] numElements;
}

void ImageUtility2::L2Check(int x, int y, int dx, int dy, Image2<int>& xNear,
    Image2<int>& yNear, Image2<int>& dist)
{
    int const dim0 = dist.GetDimension(0);
    int const dim1 = dist.GetDimension(1);
    int xp = x + dx, yp = y + dy;
    if (0 <= xp && xp < dim0 && 0 <= yp && yp < dim1)
    {
        if (dist(xp, yp) < dist(x, y))
        {
            int dx0 = xNear(xp, yp) - x;
            int dy0 = yNear(xp, yp) - y;
            int newDist = dx0*dx0 + dy0*dy0;
            if (newDist < dist(x, y))
            {
                xNear(x, y) = xNear(xp, yp);
                yNear(x, y) = yNear(xp, yp);
                dist(x, y) = newDist;
            }
        }
    }
}

bool ImageUtility2::Interior2(Image2<int>& image, int x, int y)
{
    bool b1 = (image(x, y-1) != 0);
    bool b3 = (image(x+1, y) != 0);
    bool b5 = (image(x, y+1) != 0);
    bool b7 = (image(x-1, y) != 0);
    return (b1 && b3) || (b3 && b5) || (b5 && b7) || (b7 && b1);
}

bool ImageUtility2::Interior3(Image2<int>& image, int x, int y)
{
    int numNeighbors = 0;
    if (image(x-1, y) != 0)
    {
        ++numNeighbors;
    }
    if (image(x+1, y) != 0)
    {
        ++numNeighbors;
    }
    if (image(x, y-1) != 0)
    {
        ++numNeighbors;
    }
    if (image(x, y+1) != 0)
    {
        ++numNeighbors;
    }
    return numNeighbors == 3;
}

bool ImageUtility2::Interior4(Image2<int>& image, int x, int y)
{
    return image(x-1, y) != 0
        && image(x+1, y) != 0
        && image(x, y-1) != 0
        && image(x, y+1) != 0;
}

bool ImageUtility2::MarkInterior(Image2<int>& image, int value,
    bool (*function)(Image2<int>&,int,int))
{
    int const dim0 = image.GetDimension(0);
    int const dim1 = image.GetDimension(1);
    bool noInterior = true;
    for (int y = 0; y < dim1; ++y)
    {
        for (int x = 0; x < dim0; ++x)
        {
            if (image(x, y) > 0)
            {
                if (function(image, x, y))
                {
                    image(x, y) = value;
                    noInterior = false;
                }
                else
                {
                    image(x, y) = 1;
                }
            }
        }
    }
    return noInterior;
}

bool ImageUtility2::IsArticulation(Image2<int>& image, int x, int y)
{
    // Converts 8 neighbors of pixel (x,y) to an 8-bit value, bit = 1 iff
    // pixel is set.
    int byteMask = 0;
    if (image(x-1, y-1) != 0)
    {
        byteMask |= 0x01;
    }
    if (image(x, y-1) != 0)
    {
        byteMask |= 0x02;
    }
    if (image(x+1, y-1) != 0)
    {
        byteMask |= 0x04;
    }
    if (image(x+1, y) != 0)
    {
        byteMask |= 0x08;
    }
    if (image(x+1, y+1) != 0)
    {
        byteMask |= 0x10;
    }
    if (image(x, y+1) != 0)
    {
        byteMask |= 0x20;
    }
    if (image(x-1, y+1) != 0)
    {
        byteMask |= 0x40;
    }
    if (image(x-1, y) != 0)
    {
        byteMask |= 0x80;
    }

    return msArticulation[byteMask] == 1;
}

bool ImageUtility2::ClearInteriorAdjacent(Image2<int>& image, int value)
{
    int const dim0 = image.GetDimension(0);
    int const dim1 = image.GetDimension(1);
    bool noRemoval = true;
    for (int y = 0; y < dim1; ++y)
    {
        for (int x = 0; x < dim0; ++x)
        {
            if (image(x, y) == 1)
            {
                bool interiorAdjacent =
                    image(x-1, y-1) == value ||
                    image(x  , y-1) == value ||
                    image(x+1, y-1) == value ||
                    image(x+1, y  ) == value ||
                    image(x+1, y+1) == value ||
                    image(x  , y+1) == value ||
                    image(x-1, y+1) == value ||
                    image(x-1, y  ) == value;

                if (interiorAdjacent && !IsArticulation(image, x, y))
                {
                    image(x, y) = 0;
                    noRemoval = false;
                }
            }
        }
    }
    return noRemoval;
}


int const ImageUtility2::msArticulation[256] =
{
    0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,
    0,1,1,1,1,1,1,1,0,1,0,0,0,1,0,0,
    0,1,1,1,1,1,1,1,0,1,0,0,0,1,0,0,
    0,1,1,1,1,1,1,1,0,1,0,0,0,1,0,0,
    0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
    0,1,1,1,1,1,1,1,0,1,0,0,0,1,0,0,
    0,1,1,1,1,1,1,1,0,1,0,0,0,1,0,0,
    0,0,0,0,1,1,0,0,1,1,0,0,1,1,0,0,
    1,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,
    0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,1,1,0,0,1,1,0,0,1,1,0,0,
    1,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,
    0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0
};
