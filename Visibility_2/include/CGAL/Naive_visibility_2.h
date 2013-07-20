// Copyright (c) 2013 Technical University Braunschweig (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s):  Kan Huang <huangkandiy@gmail.com>
//             
#ifndef CGAL_NAIVE_VISIBILITY_2_H
#define CGAL_NAIVE_VISIBILITY_2_H

#include <iostream>
#include <vector>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Direction_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Ray_2.h>
#include <CGAL/tags.h>
#include <CGAL/enum.h>

namespace CGAL {

namespace Visibility_2 {

//debug
template<typename Point_handle>
void print(std::vector<Point_handle> ps){
    for (int i=0; i<ps.size(); i++)
    {
        std::cout<<ps[i]->point().x()<<","<<ps[i]->point().y()<<std::endl;
    }
}
//template <typename Point_2>
//void print_point(const Point_2& p) {
//    std::cout<<"["<<p.x()<<","<<p.y()<<"]"<<std::endl;
//}


template <typename Arrangement_2, typename Regularization_tag>
class Naive_visibility_2 {
    typedef typename Arrangement_2::Geometry_traits_2         Geometry_traits_2;
    typedef typename Geometry_traits_2::Point_2						Point_2;
    typedef typename Geometry_traits_2::Ray_2						Ray_2;
    typedef typename Geometry_traits_2::Segment_2					Segment_2;
    typedef typename Geometry_traits_2::Vector_2                    Vector_2;
    typedef typename Geometry_traits_2::Direction_2                 Direction_2;

    typedef typename Arrangement_2::Halfedge             Halfedge;
    typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Arrangement_2::Vertex_const_handle  Vertex_const_handle;
    typedef typename Arrangement_2::Face_const_handle    Face_const_handle;

    enum Intersection_type { UNBOUNDED, CORNER, INNER };

public:
    //members
    Arrangement_2   arr;

    //functions
    Naive_visibility_2(const Arrangement_2 &arr):arr(arr), attach_tag(true) {}

    Naive_visibility_2(): attach_tag(false) {}


    void visibility_region(const Point_2 &q, const Halfedge &edge, Arrangement_2 &out_arr) {

    }

    void visibility_region(const Point_2 &query, Face_const_handle fh, Arrangement_2 &out_arr) {
        visibility_region_impl(query, fh, out_arr, Regularization_tag());

    }

    void visibility_region_impl(const Point_2 &query, Face_const_handle fh, Arrangement_2 &out_arr, CGAL::Tag_true) {
        std::vector<Vertex_const_handle> vertices;                    //all vertices of the face.
        std::vector<Halfedge_const_handle> edges, active_edges;       //edges stores all halfedges of the face; and active_edges stores all halfedges that is currently intersected by the view ray.
        //preprocess the face
        input_face(fh, vertices, edges, query);
        //initiation of vision ray
        Vector_2 dir;
        if (Direction_2(-1, 0) < Direction_2(Vector_2(query, (*vertices.rbegin())->point())))
        {
            dir = Vector_2(1, 0) + Vector_2(query, (*vertices.rbegin())->point());
        }
        else {
            dir = Vector_2(0, -1);
        }
        Ray_2 init_vision_ray(query, dir);

        //initiation of active_edges
        typename std::vector<Halfedge_const_handle>::iterator iter1;
        for (iter1 = edges.begin(); iter1 != edges.end(); iter1++)
        {
            insert_halfedge(active_edges, init_vision_ray, *iter1);
        }

        //angular sweep
        std::vector<Point_2> polygon;
        Ray_2 curr_vision_ray = init_vision_ray;
        typename std::vector<Vertex_const_handle>::iterator vit = vertices.begin(), begin_it, end_it;
        Halfedge_const_handle closest_edge;
        bool is_init_empty = active_edges.empty();
        while (vit != vertices.end())
        {
            if (active_edges.empty())
            {
                begin_it = vit;
                end_it = vit;
                curr_vision_ray= Ray_2(query, (*begin_it)->point());
                Direction_2 d1(curr_vision_ray), d2;
                do {
                    d2 = Direction_2(Ray_2(query, (*end_it)->point()));
                    if (d1 != d2) break;
                } while (++end_it != vertices.end());
                add_edges(begin_it, end_it, active_edges, curr_vision_ray);

                //since active_edges is empty, that means adding new edges will bring in an unbounded edge to arrangement.
                closest_edge = active_edges[0];
                remove_edges(active_edges, curr_vision_ray);

                if (active_edges.empty()) {
                    //this means in input, there is no edge that intersects curr_vision_ray
                    //except those collinear ones.
                    //because one unbounded ray has been inserted before,
                    //we don't need to do anything here.
                }
                else {
                    //add new edge;
                    Point_2 p1 = intersection_point(curr_vision_ray, closest_edge);
                    Point_2 p2 = intersection_point(curr_vision_ray, active_edges[0]);
                    update_visibility(p1, polygon, out_arr);
                    update_visibility(p2, polygon, out_arr);

                }
            }
            else {
                closest_edge = active_edges[0];
                begin_it = vit;  //all vertices between begin_it and end_it(not included) are collinear with query point.
                end_it = vit;
                curr_vision_ray = Ray_2(query, (*begin_it)->point());
                Direction_2 d1(curr_vision_ray), d2;
                do {
                    d2 = Direction_2(Ray_2(query, (*end_it)->point()));
                    if (d1 != d2) break;
                } while (++end_it != vertices.end());
                add_edges(begin_it, end_it, active_edges, curr_vision_ray);

                if (closest_edge != active_edges[0])
                {
                    //add new edge;
                    Point_2 p1 = intersection_point(curr_vision_ray, closest_edge);
                    Point_2 p2 = intersection_point(curr_vision_ray, active_edges[0]);
                    update_visibility(p1, polygon, out_arr);
                    update_visibility(p2, polygon, out_arr);

                }

                closest_edge = active_edges[0];
                remove_edges(active_edges, curr_vision_ray);

                if (active_edges.empty()) {
                    //add an unbounded edge
//todo                    CGAL::insert_curve(out_arr, Ray_2((*begin_it)->point(), Direction_2(curr_vision_ray)));
                }
                else {
                    if (closest_edge != active_edges[0]) {
                        //add new edge;
                        Point_2 p1 = intersection_point(curr_vision_ray, closest_edge);
                        Point_2 p2 = intersection_point(curr_vision_ray, active_edges[0]);
                        update_visibility(p1, polygon, out_arr);
                        update_visibility(p2, polygon, out_arr);
                    }
                }
            }
            vit = end_it;
        }
        if (!is_init_empty) {
            CGAL::insert(out_arr, Segment_2(polygon[0], polygon.back()));
        }
    }

    void visibility_region_impl(const Point_2 &query, Face_const_handle fh, Arrangement_2 &out_arr, CGAL::Tag_false) {
        std::vector<Vertex_const_handle> vertices;                    //all vertices of the face.
        std::vector<Halfedge_const_handle> edges, active_edges;       //edges stores all halfedges of the face; and active_edges stores all halfedges that is currently intersected by the view ray.
        //preprocess the face
        input_face(fh, vertices, edges, query);
        //initiation of vision ray
        Vector_2 dir;
        if (Direction_2(-1, 0) < Direction_2(Vector_2(query, (*vertices.rbegin())->point())))
        {
            dir = Vector_2(1, 0) + Vector_2(query, (*vertices.rbegin())->point());
        }
        else {
            dir = Vector_2(0, -1);
        }
        Ray_2 init_vision_ray(query, dir);

        //initiation of active_edges
        typename std::vector<Halfedge_const_handle>::iterator iter1;
        for (iter1 = edges.begin(); iter1 != edges.end(); iter1++)
        {
            insert_halfedge(active_edges, init_vision_ray, *iter1);
        }

        //angular sweep
        std::vector<Point_2> polygon;
        Ray_2 curr_vision_ray = init_vision_ray;
        typename std::vector<Vertex_const_handle>::iterator vit = vertices.begin(), begin_it, end_it;
        Halfedge_const_handle closest_edge;
        bool is_init_empty = active_edges.empty();
        while (vit != vertices.end())
        {
            if (active_edges.empty())
            {
                begin_it = vit;
                end_it = vit;
                curr_vision_ray= Ray_2(query, (*begin_it)->point());
                Direction_2 d1(curr_vision_ray), d2;
                do {
                    d2 = Direction_2(Ray_2(query, (*end_it)->point()));
                    if (d1 != d2) break;
                } while (++end_it != vertices.end());
                add_edges(begin_it, end_it, active_edges, curr_vision_ray);

                //since active_edges is empty, that means adding new edges will bring in an unbounded edge to arrangement.
                Point_2 p1 = intersection_point(curr_vision_ray, active_edges[0]);
                polygon.push_back(p1);
//todo                CGAL::insert_curve(out_arr, Ray_2(p1, d1));

                closest_edge = active_edges[0];
                remove_edges(active_edges, curr_vision_ray);

                if (active_edges.empty()) {
                    //this means there is no edge that intersects curr_vision_ray
                    //except those collinear ones.
                    //because one unbounded ray has been inserted before,
                    //we don't need to do anything here.
                }
                else {
                    //add new edge;
                    //Point_2 p1 = intersection_point(curr_vision_ray, closest_edge);
                    Point_2 p2 = intersection_point(curr_vision_ray, active_edges[0]);
                    //update_visibility(p1, polygon, out_arr);
                    update_visibility(p2, polygon, out_arr);

                }
            }
            else {
                Point_2 right_p, left_p, mid_p;
                begin_it = vit;
                end_it = vit;  //all vertices between begin_it and end_it(not included) are collinear with query point.
                curr_vision_ray= Ray_2(query, (*begin_it)->point());
                right_p = intersection_point(curr_vision_ray, active_edges[0]);
                Direction_2 d1(curr_vision_ray), d2;
                do {
                    d2 = Direction_2(Ray_2(query, (*end_it)->point()));
                    if (d1 != d2) break;
                } while (++end_it != vertices.end());
                add_edges(begin_it, end_it, active_edges, curr_vision_ray);

                mid_p = intersection_point(curr_vision_ray, active_edges[0]);
                std::vector<Point_2> collinear_vertices;
                Intersection_type i_type = needle(active_edges, curr_vision_ray, collinear_vertices);
                switch (i_type) {
                case UNBOUNDED :
                    //todo:this part is not finished.
                    //remove right and collinear;
                    remove_edges(active_edges, curr_vision_ray);
                    update_visibility(right_p, polygon, out_arr);
                    update_visibility(mid_p, polygon, out_arr);
                    //todo CGAL::insert_curve();
                    if (!active_edges.empty()) {
                        left_p = intersection_point(curr_vision_ray, active_edges[0]);
                        update_visibility(left_p, polygon, out_arr);
                    }
                    break;
                case CORNER :
                    //remove right and collinear;
                    remove_edges(active_edges, curr_vision_ray);
                    left_p = intersection_point(curr_vision_ray, active_edges[0]);
                    update_visibility(right_p, polygon, out_arr);
                    insert_needle(collinear_vertices, polygon, out_arr);
//                    update_visibility(mid_p, polygon, out_arr);
//                    update_visibility(left_p, polygon, out_arr);
                    polygon.push_back(left_p);
                    break;
                case INNER :
                    //remove right and collinear;
                    remove_edges(active_edges, curr_vision_ray);
                    if (collinear_vertices.size() < 2) {
                        //this means mid_p = left_p = right_p = furthest_p. no new vertex is found.
                    }
                    else {
                        left_p = intersection_point(curr_vision_ray, active_edges[0]);
                        update_visibility(right_p, polygon, out_arr);
                        insert_needle(collinear_vertices, polygon, out_arr);
//                        update_visibility(left_p, polygon, out_arr);
                        polygon.push_back(left_p);
                    }
                    break;
                }
            }
            vit = end_it;
        }
        if (!is_init_empty) {
            CGAL::insert(out_arr, Segment_2(polygon.back(),polygon.front()));
        }
    }


    bool is_attached() {
        return attach_tag;
    }

    void attach(Arrangement_2 arr) {
        this->arr = arr;
        this->attach_tag = true;
    }

    void detach() {
        attach_tag = false;
    }


private:

    bool            attach_tag;

    // return the intersection of a ray and a segment. if the intersection is a segment, return the end closer to the source of ray.
    // if there is no intersection, return the source of ray.
    Point_2 intersection_point(Ray_2 ray, Segment_2 seg )
    {

        CGAL::Object result = CGAL::intersection(ray, seg);
        if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
            return *ipoint;
        } else
            //if result is a segment, return the end closer to the source of ray.
            if (const Segment_2 *iseg = CGAL::object_cast<Segment_2 >(&result)) {
                switch (CGAL::compare_distance_to_point(ray.source(), iseg->source(), iseg->target())) {
                case (CGAL::SMALLER):
                    return iseg->source();
                    break;
                case (CGAL::LARGER) :
                    return iseg->target();
                    break;
                }

            } else {
                // if no intersection, return the source of ray.
                return ray.source();
            }
    }

    Point_2 intersection_point(Ray_2 ray, Halfedge_const_handle seg) {
        return intersection_point(ray, halfedge2seg(seg));
    }

    //convertor for halfedge to segment
    Segment_2 halfedge2seg(Halfedge_const_handle e){
        return Segment_2(e->source()->point(), e->target()->point());
    }

    //check whether two halfedges are the same segment.
    bool is_same_edge(Halfedge_const_handle e1, Halfedge_const_handle e2) {

    }

    //given two edges incident to a vision ray at the same point, find which one is first seen in sweeping.
    bool is_closer(const Ray_2 &ray, Halfedge_const_handle seg1, Halfedge_const_handle seg2) {
        Point_2 shared = intersection_point(ray, seg1);
        Point_2 end1, end2;
        if (shared == seg1->source()->point())
            end1 = seg1->target()->point();
        else
            end1 = seg1->source()->point();

        if (shared == seg2->source()->point())
            end2 = seg2->target()->point();
        else
            end2 = seg2->source()->point();
        if (CGAL::right_turn(ray.source(), shared, end1) && !CGAL::right_turn(ray.source(), shared, end2))
            return true;
        if (CGAL::right_turn(ray.source(), shared, end2) && !CGAL::right_turn(ray.source(), shared, end1))
            return false;
        switch (CGAL::orientation(ray.source(), shared, end1)) {
        case CGAL::COLLINEAR:
            return (CGAL::left_turn(ray.source(), shared, end2));
        case CGAL::RIGHT_TURN:
            return (CGAL::right_turn(end1, shared, end2));
        case CGAL::LEFT_TURN:
            return (CGAL::left_turn(end1, shared, end2));
        }
    }

    //insert newly-discovered edges into active_edges according to its intersection with the view ray.
    void insert_halfedge(std::vector<Halfedge_const_handle> &active_edges, const Ray_2 &ray, Halfedge_const_handle edge)
    {
        Point_2 cross_of_e = intersection_point(ray, edge);
        if (cross_of_e != ray.source())
        {
            typename std::vector<Halfedge_const_handle>::iterator curr = active_edges.begin();
            while (curr != active_edges.end())
            {

                Point_2 cross_of_curr = intersection_point(ray, *curr);
                if (CGAL::compare_distance_to_point(ray.source(), cross_of_e, cross_of_curr) == CGAL::SMALLER)
                    break;                
                if (cross_of_curr == cross_of_e && is_closer(ray, edge, *curr))
                    break;
                ++curr;
            }
            active_edges.insert(curr, edge);
        }
    }

    //insert vh into vertices by the angle between ray<p, vh> and positive x-ray.
    //if the angles are the same, compare their distances to p.
    void sort_vertex(std::vector<Vertex_const_handle>& vertices, Vertex_const_handle vh, const Point_2& p)
    {
        typename std::vector<Vertex_const_handle>::iterator first = vertices.begin();
        Vector_2 vector_of_v(p, vh->point());
        Direction_2 dir_of_v(vector_of_v);
        while (first != vertices.end())
        {
            if (vh->point() == (*first)->point()) {
                return;
            }
            Vector_2 vector1(p, (*first)->point());
            Direction_2 d1(vector1);
            if (dir_of_v < d1)
                break;
            //if same angles are the same, put the vertex which is closer to query point before.
            if (dir_of_v == d1 && CGAL::compare_distance_to_point(p, vh->point(), (*first)->point()) == CGAL::SMALLER)
                break;
            ++first;
        }
        vertices.insert(first, vh);
    }


    //sort vertex vh by the angle between vector<p, v> and positive x-ray.
    //if the angles are the same, place them in a 'circular' way.
    //we do this because the output is not regularized, i.e. 1-d needle is reserved in output.
//    void sort_vertex(vector<Vertex_const_handle>& vertices, Vertex_const_handle vh, const Point_2& p, CGAL::Tag_false)
//    {
//        Ray_2 ray(p, vh->point());
//        typename vector<Vertex_const_handle>::iterator first = vertices.begin();
//        while (first != vertices.end())
//        {
//            Ray_2 r(p, (*first)->point());
//            Direction_2 d1(r);
//            if (ray.direction() < d1)
//                break;
//            //if angles are the same, then using Isblock() to decide the order of vertices on the view ray.
//            if (ray.direction() == d1)
//            {
//                if (Is_block(vh, ray) && (!Is_block(*first, ray)))
//                    break;
//                if (Is_block(vh, ray) && Is_block(*first, ray) && CGAL::compare_distance_to_point(p, vh->point(), (*first)->point()) == CGAL::SMALLER)
//                    break;
//                if (!Is_block(vh, ray) && !Is_block(*first, ray) && CGAL::compare_distance_to_point(p, vh->point(), (*first)->point()) == CGAL::LARGER)
//                    break;
//            }
//            ++first;
//        }
//        vertices.insert(first, vh);
//    }

    //traverse the face to get all edges and sort vertices in counter-clockwise order.
    void input_face (Face_const_handle fh,
                     std::vector<Vertex_const_handle>& vertices,
                     std::vector<Halfedge_const_handle>& edges,
                     const Point_2& p)
    {
        typename Arrangement_2::Ccb_halfedge_const_circulator curr = fh->outer_ccb();
        typename Arrangement_2::Ccb_halfedge_const_circulator circ = curr;
        do {
            sort_vertex(vertices, curr->source(), p);
            edges.push_back(curr);
        } while (++curr != circ);
        typename Arrangement_2::Hole_const_iterator hi;
        for (hi = fh->holes_begin(); hi != fh->holes_end(); ++hi) {
            typename Arrangement_2::Ccb_halfedge_const_circulator c1 = *hi, c2 = *hi;
            do {
                sort_vertex(vertices, c1->source(), p);
                edges.push_back(c1);
            } while (++c1 != c2);
        }
    }




    //insert new vertice to polygon. before insertion, check if this vertice has been added before.
    void update_visibility(const Point_2 p, std::vector<Point_2>& polygon, Arrangement_2& arr){
        if (polygon.empty())
            polygon.push_back(p);
        else
        {
            if (polygon.back() != p){
                CGAL::insert(arr, Segment_2(polygon.back(), p));
                polygon.push_back(p);
            }
        }
    }

    void insert_needle(const std::vector<Point_2>& points, std::vector<Point_2>& polygon, Arrangement_2 &arr){
        if (points.size() > 1) {
            for (int i = 0; i != points.size()-1; i++) {
                CGAL::insert(arr, Segment_2(points[i], points[i+1]));
            }
        }
    }


    //add a new edge when vision ray passes a vertex
    void add_edge(Vertex_const_handle vh,
                   std::vector<Halfedge_const_handle>& edges,
                   const Ray_2& r) {
        typename Arrangement_2::Halfedge_around_vertex_const_circulator first, curr;
        first = curr = vh->incident_halfedges();
        do {
            if (!CGAL::right_turn(r.source(), vh->point(), curr->source()->point()))
            {
                insert_halfedge(edges, r, curr);
            }
        } while (++curr != first);

    }
    //add new edges
    void add_edges(typename std::vector<Vertex_const_handle>::iterator begin_it,
                   typename std::vector<Vertex_const_handle>::iterator end_it,
                   std::vector<Halfedge_const_handle>& edges,
                   const Ray_2& r)
    {
        do {
            add_edge(*begin_it, edges, r);
        } while (++begin_it != end_it);

    }

    //remove edges that are not active any longer
    void remove_edges(std::vector<Halfedge_const_handle>& edges, const Ray_2& r) {
        typename std::vector<Halfedge_const_handle>::iterator eit = edges.begin();
        while (eit != edges.end()) {
            Point_2 p1 = (*eit)->target()->point();
            Point_2 p2 = (*eit)->source()->point();
            bool is_incident(false);
            if (is_on_ray(r, p1)) {
                is_incident = true;
            }
            else if (is_on_ray(r, p2)) {
                Point_2 tmp = p1;
                p1 = p2;
                p2 = p1;
                is_incident = true;
            }


            if ( is_incident && !CGAL::left_turn(r.source(), p1, p2))
            {
                eit = edges.erase(eit);
            }
            else {
                eit++;
            }
        }


    }

    bool is_on_ray(const Ray_2& r, const Point_2& p) {
        return Direction_2(Vector_2(r.source(), p)) == Direction_2(r);
    }
    //return the type of the needle.
    //the vertices on the needle will be saved in collinear_vertices.
    Intersection_type needle(std::vector<Halfedge_const_handle>& edges, Ray_2& r, std::vector<Point_2>& collinear_vertices) {
        typename std::vector<Halfedge_const_handle>::iterator curr = edges.begin();
//        Point_2 p = r.source(), end1, end2;
        Vertex_const_handle vertex1;
        //flag shows whether the left side or right side of needle is blocked.
        bool block_left, block_right;
        do {
            Point_2 cross = intersection_point(r, *curr);
            if (cross != (*curr)->source()->point() && cross != (*curr)->target()->point()) {
                collinear_vertices.push_back(cross);
                return INNER;
            }
            if (CGAL::orientation(r.source(), (*curr)->source()->point(), (*curr)->target()->point()) == CGAL::COLLINEAR) {
                vertex1 = (*curr)->target();
            }
            else {
                if (cross == (*curr)->source()->point()) {
                    vertex1 = (*curr)->source();
                }
                else {
                    vertex1 = (*curr)->target();
                }
            }
            if (collinear_vertices.empty() || vertex1->point() != collinear_vertices.back()) {
                collinear_vertices.push_back(vertex1->point());
                //flag shows whether the left side or right side of current vertex is blocked.
                //has_predecessor indicates whether this vertex is incident to an edge whose another end is between the source of ray and it,
                //because that will effect the value of block_left, block_right.
                bool left_v(false), right_v(false), has_predecessor(false);

                typename Arrangement_2::Halfedge_around_vertex_const_circulator first_edge, curr_edge;
                first_edge = curr_edge = vertex1->incident_halfedges();
                do {
                    switch (CGAL::orientation(r.source(), curr_edge->target()->point(), curr_edge->source()->point())) {
                    case CGAL::RIGHT_TURN :
                        right_v = true;
                        break;
                    case CGAL::LEFT_TURN :
                        left_v = true;
                        break;
                    case CGAL::COLLINEAR :
                        if (CGAL::compare_distance_to_point(r.source(), curr_edge->target()->point(), curr_edge->source()->point()) == CGAL::LARGER) {
                            has_predecessor = true;
                        }
                    }

                } while (++curr_edge != first_edge);
                if (has_predecessor) {
                    block_left = block_left || left_v;
                    block_right = block_right || right_v;
                }
                else {
                    block_left = left_v;
                    block_right = right_v;
                }
                if (block_left && block_right) {
                    return CORNER;
                }
            }
        } while (++curr != edges.end());
        return UNBOUNDED;
    }
//debug
    void print_edges(std::vector<Halfedge_const_handle>& edges){
        for (int i = 0; i != edges.size(); i++) {
            Point_2 p1, p2;
            p1 = edges[i]->source()->point();
            p2 = edges[i]->target()->point();
            std::cout<<p1<<"->"<<p2<<std::endl;
        }
    }


    //angular sweep a vertice of face.
    void sweep_vertex(std::vector<Halfedge_const_handle> &active_edges, const Point_2 &query, Vertex_const_handle vh, std::vector<Point_2> &polygon )
    {
        //closest_edge_copy is a copy of the closest edge to query point in active_edges before sweeping.
        Halfedge_const_handle closest_edge_copy = active_edges[0];
        Ray_2 ray(query, vh->point());
        int add_count(0);
        int del_count(0);

        //delete all edges in active_edges which is incident to v, because they has been sweeped over
        typename std::vector<Halfedge_const_handle>::iterator edge_iter = active_edges.begin();
        while (edge_iter != active_edges.end()) {
            if (((*edge_iter)->source()->point() == vh->point()) || ((*edge_iter)->target()->point() == vh->point()))
            {
                edge_iter = active_edges.erase(edge_iter);
                ++del_count;
            }
            else {
                ++edge_iter;
            }
        }

        //add all edges which is incident to v but not in active_edges before to active_edges
        typename Arrangement_2::Halfedge_around_vertex_const_circulator first, curr;
        first = curr = vh->incident_halfedges();
        do {
            if (CGAL::left_turn(query, vh->point(), curr->source()->point()))
            {
                insert_halfedge(active_edges, ray, curr);
                ++add_count;
            }
            else if (CGAL::collinear(query, vh->point(), curr->source()->point()) &&
                    CGAL::compare_distance_to_point(query, vh->point(), curr->source()->point()) == CGAL::SMALLER)
            {
                insert_halfedge(active_edges, ray, curr);
                ++add_count;
            }
        } while (++curr != first);

        //update the visibility region
        if (closest_edge_copy != active_edges[0]) {
            //when the closest edge changed
            if (del_count > 0 && add_count > 0) {
                //some edges are added and some are deleted, which means the vertice sweeped is a vertice of visibility polygon.
                update_visibility(vh->point(), polygon);
            }
            if (del_count == 0 && add_count > 0) {
                //only add some edges, means the view ray is blocked by new edges.
                //therefore first add the intersection of view ray and former closet edge, then add the vertice sweeped.
                update_visibility(intersection_point(ray, halfedge2seg(closest_edge_copy)), polygon);
                update_visibility(vh->point(), polygon);
            }
            if (del_count > 0 && add_count == 0) {
                //only delete some edges, means some block is moved and the view ray can reach the segments after the block.
                update_visibility(vh->point(), polygon);
                update_visibility(intersection_point(ray, halfedge2seg(active_edges[0])), polygon);
            }
        }

    }

};

//For debug. Print all edges of arrangements into console.
template <typename Arrangement_2>
void print_arrangement(const Arrangement_2 &arr) {
    typedef typename Arrangement_2::Edge_const_iterator Edge_const_iterator;
    Edge_const_iterator eit;
    std::cout << arr.number_of_edges() << " edges:" << std::endl;
    for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit)
      std::cout << "[" << eit->curve() << "]" << std::endl;
}



} // namespace Visibility_2
} // namespace CGAL



#endif


