#ifndef BASICANNKDTREE_H
#define BASICANNKDTREE_H

#include <ANN/ANN.h>
#include <vector>

class BasicANNkdTree
{
    unsigned int nb_points;
    ANNpointArray ANNvertex;
    ANNkd_tree *ANNtree;
    size_t points_dimension;

public:
    BasicANNkdTree() : nb_points( 0 ) , ANNvertex( NULL ) , ANNtree( NULL ) , points_dimension(3) {}

    ANNkd_tree * tree() { return ANNtree; }

    void clear() {
        if( ANNtree != NULL ) {
            delete ANNtree;
            ANNtree = NULL;
        }
        if( ANNvertex != NULL ) {
            annDeallocPts(ANNvertex);
            ANNvertex = NULL;
        }
    }

    template< class point_t >
    void build( std::vector< point_t > const & i_ps ) {
        if(ANNvertex != NULL) {
            annDeallocPts(ANNvertex);
            ANNvertex = NULL;
        }
        ANNvertex = annAllocPts(i_ps.size(),points_dimension);

        for(unsigned int v = 0 ; v < i_ps.size() ; ++v ) {
            point_t const & p = i_ps[v];
            for(unsigned int dimIt = 0 ; dimIt < points_dimension ; ++dimIt ) {
                ANNvertex[v][dimIt] = p[dimIt];
            }
        }
        if(ANNtree != NULL)
            delete ANNtree;

        ANNtree = new ANNkd_tree(ANNvertex,i_ps.size(),points_dimension,1);
    }






    template< class point_t >
    inline unsigned int nearest( point_t const & i_position ) const {
        ANNpoint ann_point = annAllocPt(points_dimension);
        for(unsigned int dimIt = 0 ; dimIt < points_dimension ; ++dimIt ) {
            ann_point[dimIt] = i_position[dimIt];
        }
        ANNidx idx ; ANNdist dd;
        ANNtree->annkSearch( ann_point , 1 , &idx , &dd );
        annDeallocPt(ann_point);
        return (unsigned int)( idx );
    }




    // BEFORE you call the following function (knearest), be careful , you need to allocate first id_nearest_neighbors and square_distances_to_neighbors this way :
    //      ANNidxArray id_nearest_neighbors = new ANNidx[ k ]
    //      ANNdistArray square_distances_to_neighbors = new ANNdist[ k ]
    // Don't forget after to delete your data when you're done :
    //      delete [] id_nearest_neighbors;
    //      delete [] square_distances_to_neighbors;
    // It gives you in id_voisins the id of the k nearest point in the tree to i_position, they are given by
    // increasing distance to the input point i_position.  ( voisins_square_distances is the vector of the SQUARE distances !! )
    template< class point_t >
    inline
            void knearest( point_t const & i_position , int k , ANNidxArray id_nearest_neighbors , ANNdistArray square_distances_to_neighbors ) const {
        ANNpoint ann_point = annAllocPt(points_dimension);
        for(unsigned int dimIt = 0 ; dimIt < points_dimension ; ++dimIt ) {
            ann_point[dimIt] = i_position[dimIt];
        }
        ANNtree->annkSearch( ann_point , k , id_nearest_neighbors , square_distances_to_neighbors );
        annDeallocPt(ann_point);
    }
};


#endif // BASICANNKDTREE_H
