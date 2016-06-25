// vim: tabstop=4:softtabstop=4:shiftwidth=4:expandtab:smarttab
#import <Foundation/Foundation.h>

@interface NNetwork : NSObject
{
@private
    double x;
    double y;
}

- (id) x: (double) x_value;
- (double) x;
- (id) y: (double) y_value;
- (double) y;
- (double) magnitude;
@end

