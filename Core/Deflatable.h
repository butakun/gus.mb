// $Id: Deflatable.h 76 2010-10-06 16:50:49Z kato $
#ifndef INCLUDED_DEFLATABLE_H__
#define INCLUDED_DEFLATABLE_H__

#include <cstdlib>

class Deflatable
{
public:
    /// Returns the byte size of the deflated object, if deflated.
    virtual size_t DeflatedSize() const = 0;

    /// Deflates the object, and returns the byte size of the deflated object
    virtual size_t Deflate(char* data, size_t buflen) = 0;

    /// Inflates the data into an object
    virtual void Inflate(const char* data, size_t buflen) = 0;

protected:

private:
};

#endif // INCLUDED_DEFLATABLE_H__

