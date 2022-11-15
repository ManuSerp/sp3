#include <assert.h>
#include <stddef.h>

#include "heatsim.h"
#include "log.h"

struct parametres_s {
    unsigned int width;
    unsigned int height;
    int padding;

};

struct data_s {

    double *data;

};

static MPI_Datatype parametres_type;
static MPI_Datatype data_type;




int heatsim_init(heatsim_t* heatsim, unsigned int dim_x, unsigned int dim_y) {
    /*
     * TODO: Initialiser tous les membres de la structure `heatsim`.
     *       Le communicateur doit être périodique. Le communicateur
     *       cartésien est périodique en X et Y.
     */

    int dims[2] = {dim_x, dim_y}; // peut etre unsigned int
    int periods[2] = {1, 0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1,&heatsim->communicator);

    MPI_Comm_rank(heatsim->communicator, &heatsim->rank);
    MPI_Comm_size(heatsim->communicator,&heatsim->rank_count);

    MPI_Cart_coords(heatsim->communicator, heatsim->rank, 2, heatsim->coordinates);

    int north[2]={heatsim->coordinates[0], heatsim->coordinates[1] - 1};
    int south[2]={heatsim->coordinates[0], heatsim->coordinates[1] + 1};
    int east[2]={heatsim->coordinates[0] + 1, heatsim->coordinates[1]};
    int west[2]={heatsim->coordinates[0] - 1, heatsim->coordinates[1]};



    MPI_Cart_rank(heatsim->communicator,north, &heatsim->rank_north_peer);
    MPI_Cart_rank(heatsim->communicator,south, &heatsim->rank_south_peer);
    MPI_Cart_rank(heatsim->communicator,east, &heatsim->rank_east_peer);
    int error = MPI_Cart_rank(heatsim->communicator,west, &heatsim->rank_west_peer);
    if (error != MPI_SUCCESS) {
        goto fail_exit;   

  }

  //param datatype
    MPI_Datatype p_fields[3] = {MPI_UNSIGNED, MPI_UNSIGNED, MPI_INT};

    int p_blocklen[3] = {1, 1, 1};

    MPI_Aint p_displ[3];
    p_displ[0] = offsetof(struct parametres_s, width);
    p_displ[1] = offsetof(struct parametres_s, height);
    p_displ[2] = offsetof(struct parametres_s, padding);

    MPI_Type_create_struct(3, p_blocklen, p_displ, p_fields, &parametres_type);

    MPI_Type_commit(&parametres_type);

    // data datatype

    MPI_Datatype d_fields[1] = {MPI_DOUBLE};
    int d_blocklen[1] = {1};
    MPI_Aint d_displ[1];
    d_displ[0] = offsetof(struct data_s, data);

    MPI_Type_create_struct(1, d_blocklen, d_displ, d_fields, &data_type);

    MPI_Type_commit(&data_type);



fail_exit:
    return -1;
}

int heatsim_send_grids(heatsim_t* heatsim, cart2d_t* cart) {
    /*
     * TODO: Envoyer toutes les `grid` aux autres rangs. Cette fonction
     *       est appelé pour le rang 0. Par exemple, si le rang 3 est à la
     *       coordonnée cartésienne (0, 2), alors on y envoit le `grid`
     *       à la position (0, 2) dans `cart`.
     *
     *       Il est recommandé d'envoyer les paramètres `width`, `height`
     *       et `padding` avant les données. De cette manière, le receveur
     *       peut allouer la structure avec `grid_create` directement.
     *
     *       Utilisez `cart2d_get_grid` pour obtenir la `grid` à une coordonnée.
     */
    struct parametres_s p_t;
    struct data_s d_t;
    
    for (int i = 1; i < heatsim->rank_count; i++ ){
            int coords[2];
            MPI_Cart_coords(heatsim->communicator, i, 2, coords);
            grid_t *grid = cart2d_get_grid(cart, coords[0], coords[1]);
            
            p_t.width = grid->width;
            p_t.height = grid->height;
            p_t.padding = grid->padding;
            d_t.data = grid->data;

            MPI_Request request;

            MPI_Isend(&p_t, 1, parametres_type, i, 0, heatsim->communicator, &request);
            MPI_Wait(&request, MPI_STATUS_IGNORE);
            MPI_Isend(&d_t, 1, data_type, i, 0, heatsim->communicator, &request);

        
    }




return 0;
}

grid_t* heatsim_receive_grid(heatsim_t* heatsim) {
    /*
     * TODO: Recevoir un `grid ` du rang 0. Il est important de noté que
     *       toutes les `grid` ne sont pas nécessairement de la même
     *       dimension (habituellement ±1 en largeur et hauteur). Utilisez
     *       la fonction `grid_create` pour allouer un `grid`.
     *
     *       Utilisez `grid_create` pour allouer le `grid` à retourner.
     */

    struct parametres_s p_t;
    struct data_s d_t;

    MPI_Request request1;
    MPI_Request request2;


    MPI_Irecv(&p_t, 1, parametres_type, 0, 0, heatsim->communicator, &request1);
    MPI_Wait(&request1, MPI_STATUS_IGNORE);

    MPI_Irecv(&d_t, p_t->width*p_t->height, data_type, 0, 0, heatsim->communicator, &request2);
    MPI_Wait(&request2, MPI_STATUS_IGNORE);


    grid_t *grid = grid_create(p_t.width, p_t.height, p_t.padding);
    grid_set(grid, *d_t.data);

    return grid;
}

int heatsim_exchange_borders(heatsim_t* heatsim, grid_t* grid) {
    assert(grid->padding == 1);

    MPI_Request request[8];

    //send
    



    /*
     * TODO: Échange les bordures de `grid`, excluant le rembourrage, dans le
     *       rembourrage du voisin de ce rang. Par exemple, soit la `grid`
     *       4x4 suivante,
     *
     *                            +-------------+
     *                            | x x x x x x |
     *                            | x A B C D x |
     *                            | x E F G H x |
     *                            | x I J K L x |
     *                            | x M N O P x |
     *                            | x x x x x x |
     *                            +-------------+
     *
     *       où `x` est le rembourrage (padding = 1). Ce rang devrait envoyer
     *
     *        - la bordure [A B C D] au rang nord,
     *        - la bordure [M N O P] au rang sud,
     *        - la bordure [A E I M] au rang ouest et
     *        - la bordure [D H L P] au rang est.
     *
     *       Ce rang devrait aussi recevoir dans son rembourrage
     *
     *        - la bordure [A B C D] du rang sud,
     *        - la bordure [M N O P] du rang nord,
     *        - la bordure [A E I M] du rang est et
     *        - la bordure [D H L P] du rang ouest.
     *
     *       Après l'échange, le `grid` devrait avoir ces données dans son
     *       rembourrage provenant des voisins:
     *
     *                            +-------------+
     *                            | x m n o p x |
     *                            | d A B C D a |
     *                            | h E F G H e |
     *                            | l I J K L i |
     *                            | p M N O P m |
     *                            | x a b c d x |
     *                            +-------------+
     *
     *       Utilisez `grid_get_cell` pour obtenir un pointeur vers une cellule.
     */

    return 0;
}

int heatsim_send_result(heatsim_t* heatsim, grid_t* grid) {
    assert(grid->padding == 0);

    /*
     * TODO: Envoyer les données (`data`) du `grid` résultant au rang 0. Le
     *       `grid` n'a aucun rembourage (padding = 0);
     */

    return 0;
}

int heatsim_receive_results(heatsim_t* heatsim, cart2d_t* cart) {
    /*
     * TODO: Recevoir toutes les `grid` des autres rangs. Aucune `grid`
     *       n'a de rembourage (padding = 0).
     *
     *       Utilisez `cart2d_get_grid` pour obtenir la `grid` à une coordonnée
     *       qui va recevoir le contenue (`data`) d'un autre noeud.
     */

    return 0;
}
