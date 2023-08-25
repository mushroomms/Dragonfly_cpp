#include <gmpxx.h>
#include <assert.h>
#include <iostream>
#include <cstdio>
#include <random>
#include "include/point.h"
#include "include/curve.h"
#include "include/mpz_math.h"
#include "include/peer.h"

#include <arpa/inet.h>
#include <sys/socket.h>
#include <iostream>
#include <sys/ioctl.h>
#include <net/if.h>
#include <unistd.h>
#include <sodium.h>
#include <cstring>
#include <sstream>
#include <fstream>
#include <x86intrin.h>

#define cpucycles(cycles) cycles = __rdtsc()
#define cpucycles_reset() cpucycles_sum = 0
#define cpucycles_start() cpucycles(cpucycles_before)
#define cpucycles_stop()                                 \
  do                                                     \
  {                                                      \
    cpucycles(cpucycles_after);                          \
    cpucycles_sum += cpucycles_after - cpucycles_before; \
  } while (0)

#define cpucycles_result() cpucycles_sum

unsigned long long cpucycles_before, cpucycles_after, cpucycles_sum;

using namespace std;

std::string getMACAddress(const char* interfaceName, int fd) {
    struct ifreq ifr;

    std::strcpy(ifr.ifr_name, interfaceName);

    if (ioctl(fd, SIOCGIFHWADDR, &ifr) == -1) {
        std::cerr << "Error getting MAC address" << std::endl;
        return "";
    }

    char macAddress[18];
    std::snprintf(macAddress, sizeof(macAddress), "%.2X:%.2X:%.2X:%.2X:%.2X:%.2X",
        static_cast<unsigned char>(ifr.ifr_hwaddr.sa_data[0]),
        static_cast<unsigned char>(ifr.ifr_hwaddr.sa_data[1]),
        static_cast<unsigned char>(ifr.ifr_hwaddr.sa_data[2]),
        static_cast<unsigned char>(ifr.ifr_hwaddr.sa_data[3]),
        static_cast<unsigned char>(ifr.ifr_hwaddr.sa_data[4]),
        static_cast<unsigned char>(ifr.ifr_hwaddr.sa_data[5]));
    
    return macAddress;
}

int main(int argc, char* argv[]) {
    srand(time(NULL));

    if (sodium_init() != 0) {
        printf("panic! the library couldn't be initialized; it is not safe to use");
        return 1;
    }

    /**
    mpz_t a, b, c, d;
    mpz_init(a);
    mpz_init(b);
    mpz_init(c);
    mpz_init(d);
    //计算2的1000次方
    mpz_init_set_ui(a, 2);
    mpz_pow_ui(c, a, 1000);
    gmp_printf("c = %Zd\n", c);

    //计算12345678900987654321*98765432100123456789
    mpz_init_set_str(b, "41577311594322086732010735089046630069560194285817217577730709794461777884595", 10); //10进制
    mpz_init_set_str(c, "76884956397045344220809746629001649093037950200943055203735601445031516197751", 10);
    mpz_mul(d, b, c);
    gmp_printf("d = %Zd\n", d);

    //计算除法
    // mpz_init_set_str(c, "9", 10);
    // mpz_init_set_str(b, "2", 10);//10进制
    mpz_fdiv_q(d, c, b); //向下整除
    gmp_printf("mpz_fdiv_q:d = %Zd\n", d);

    legendre(d, b, c);
    gmp_printf("legendre:d = %Zd\n", d);

    tonelli_shanks(d, b, c);
    gmp_printf("tonelli_shanks:d = %Zd\n", d);

    EllipticCurve ec(
        "56698187605326110043627228396178346077120614539475214109386828188763884139993",
        "17577232497321838841075697789794520262950426058923084567046852300633325438902",
        "76884956397045344220809746629001649093037950200943055203735601445031516197751", 10);
    mpz_t x;

    mpz_init_set_str(x, "44636889021951429617878751535492807962563326887208848194489661399766546537096", 10);
    ec.curve_equation(d, x);
    gmp_printf("curve_equation :d = %Zd\n", d);
    mpz_clear(x);

    bool is_QR;
    is_QR = ec.is_quadratic_residue(d);
    gmp_printf("is_quadratic_residue : d = %d\n", is_QR);

    mpz_init_set_str(x, "1128969709662994286855259741431387894916776047019745036539814346559507564901", 10);
    is_QR = ec.is_quadratic_residue(x);
    gmp_printf("is_quadratic_residue :d = %d\n", is_QR);
    mpz_clear(x);

    Point p(
        "11245091364210913538031558554069068295821326089979728553744121109670093869553",
        "58227085528211006908619226448297745234381220210201339909209560974723630792484",
        10);
    bool is_on_curve = ec.is_point_on_curve(p);
    gmp_printf("is_point_on_curve :d = %d\n", is_on_curve);

    mpz_init_set_str(x, "31585362742042248681507141308670880101726446859422047342877633324024259170673", 10);
    ec.inv_mod_p(d, x);
    gmp_printf("inv_mod_p :d = %Zd\n", d);

    Point p1(
        "63463127379121579880278198109336113474044035466585342005071192132059818598209",
        "16558137935435870897706968960488825780425063076841601361802226641039654238867",
        10);
    Point p2;
    ec.ec_point_inv(p2, p1);
    gmp_printf("ec_point_inv :x = %Zd , y = %Zd \n", p2.x, p2.y);

    Point p3(
        "44798834989657919863089174096808959741372081687640931490811108075146543214623",
        " 32552174646444655060134058623257090394879394920791321656975553975238125448510",
        10);
    Point p4;
    ec.ec_point_double_add(p4, p3);
    gmp_printf("ec_point_double_add :x = %Zd , y = %Zd \n", p4.x, p4.y);

    p1.set_x("19255648800734979802927059587974852932026844473492005743471589635415939677617", 10);
    p1.set_y("5347230032267046714807334904413528386238164048132981468634421549716907951517", 10);
    p2.set_x("53966716554608051077465919234931511423879734331111662366752677474299058504440", 10);
    p2.set_y("8117559023792801597289705548012874557664944241595251311237145384718921880283", 10);
    ec.ec_point_add(p4, p1, p2);
    gmp_printf("ec_point_add :x = %Zd , y = %Zd \n", p4.x, p4.y);

    p3.set_x("72639184339425654576555087557111920487359235746629312591395173225124908193763", 10);
    p3.set_y("15224464906049663749133738389570345770876661096229951121126739762545982893453", 10);
    mpz_init_set_str(x, "30365056390190848852943231209253184370325221626750165764764158042136144461840", 10);
    ec.ec_point_scalar_mul(p4, x, p3);
    mpz_clear(x);
    gmp_printf("ec_point_scalar_mul :x = %Zd , y = %Zd \n", p4.x, p4.y);

    p3.set_x("1104832230923804541019855326938339287498267286713819743046729009762879975423", 10);
    p3.set_y("52921860843340567879358341355783097105123683025468323692090409904865326996090", 10);
    mpz_init_set_str(x, "25568832556651077009177958848123513550279603375538798880147678520829252774745", 10);
    ec.ec_point_scalar_mul(p4, x, p3);
    mpz_clear(x);
    gmp_printf("ec_point_scalar_mul :x = %Zd , y = %Zd \n", p4.x, p4.y);
    */

    char *SERVER_IP = argv[1];
    int PORT = atoi(argv[2]);
    if (argc != 3) {
        printf("Usage: %s <server_ip> <server_port>\n", argv[0]);
        return 1;
    }

    int status, client_fd, valread, bytes_sent;
    struct sockaddr_in serv_addr;
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_port = htons(PORT);

    if ((client_fd = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0) {
        printf("\n Socket creation error \n");
        return 1;
    }
    // Convert IPv4 and IPv6 addresses from text to binary form
    if (inet_pton(AF_INET, SERVER_IP, &serv_addr.sin_addr) <= 0) {
        printf("\nInvalid address/ Address not supported \n");
        return 1;
    }
    if ((status = connect(client_fd, (struct sockaddr*)&serv_addr, sizeof(serv_addr))) < 0) {
        printf("\nConnection Failed \n");
        return 1;
    }

    gmp_printf("---------------------------------------------------\n");

    std::string own_mac = getMACAddress("ens33", client_fd); // Change to your interface name
    Peer sta("abc1238", own_mac, "STA");

    // Start timer
    struct timespec begin_sending_mac, end_receiving_mac, begin_sending_scalar, end_receiving_scalar, begin_sending_token, end_receiving_token;
    struct timespec begin_dragonfly_wall, begin_dragonfly_cpu, end_dragonfly_wall, end_dragonfly_cpu;
    clock_gettime(CLOCK_REALTIME, &begin_dragonfly_wall);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin_dragonfly_cpu);

    uint8_t ap_mac[18];
    bytes_sent = send(client_fd, own_mac.c_str(), 18, 0);
    clock_gettime(CLOCK_REALTIME, &begin_sending_mac);
    std::cout << "[STA] Transmitted Data Size " << bytes_sent << " Bytes." << std::endl;
    std::cout << "[STA] STA MAC Address sent." << std::endl;

    valread = recv(client_fd, ap_mac, sizeof(ap_mac), 0);
    clock_gettime(CLOCK_REALTIME, &end_receiving_mac);
    std::cout << "\n[STA] Received Data Size " << valread << " Bytes." << std::endl;
    std::cout << "[STA] AP MAC received." << std::endl;

    double RTT_MAC = (end_receiving_mac.tv_sec - begin_sending_mac.tv_sec) + (end_receiving_mac.tv_nsec - begin_sending_mac.tv_nsec) / 1000000000.0 * 1000.0;
    
    std::cout << "\n[STA] Own MAC Address: " << own_mac << std::endl;
    std::cout << "[STA] Peer MAC Address: " << ap_mac << std::endl;

    gmp_printf("---------------------------------------------------\n");    
    gmp_printf("Starting hunting and pecking to derive PE...\n\n");
    std::string mac_string_ap(reinterpret_cast<char*>(ap_mac), sizeof(ap_mac));
    
    cpucycles_reset();
    cpucycles_start();
    sta.initiate(mac_string_ap);
    cpucycles_stop();
    unsigned int sta_initiate_cycles = cpucycles_result();

    gmp_printf("---------------------------------------------------\n");
    gmp_printf("Starting dragonfly commit exchange...\n\n");

    cpucycles_reset();
    cpucycles_start();
    sta.commit_exchange();
    cpucycles_stop();
    unsigned int sta_commit_cycles = cpucycles_result();

    // Putting Scalar and Element into stream to be sent to peer
    std::ostringstream scalar_element_stream;

    // Stream the scalar and element
    scalar_element_stream << mpz_get_str(nullptr, 10, sta.scalar) << "\n";
    scalar_element_stream << mpz_get_str(nullptr, 10, sta.element.x) << "\n";
    scalar_element_stream << mpz_get_str(nullptr, 10, sta.element.y);
    std::string scalar_complete_stream = scalar_element_stream.str();

    // Sending Scalar & Element to peer
    bytes_sent = send(client_fd, scalar_complete_stream.c_str(), 512, 0);
    clock_gettime(CLOCK_REALTIME, &begin_sending_scalar);
    std::cout << "\n[STA] Transmitted Data Size " << bytes_sent << " Bytes." << std::endl;
    std::cout << "[STA] STA Scalar & Element sent." << std::endl;

    // Receiving peer Scalar & Element
    uint8_t peer_scalarElement[512];
    valread = recv(client_fd, peer_scalarElement, sizeof(peer_scalarElement), 0);
    clock_gettime(CLOCK_REALTIME, &end_receiving_scalar);
    std::cout << "\n[STA] Received Data Size " << valread << " Bytes." << std::endl;
    std::cout << "[STA] AP Scalar & Element received." << std::endl;

    double RTT_SCALAR = (end_receiving_scalar.tv_sec - begin_sending_scalar.tv_sec) + (end_receiving_scalar.tv_nsec - begin_sending_scalar.tv_nsec) / 1000000000.0 * 1000.0;

    // Splitting received Scalar and Element from stream
    mpz_t scalar_ap, element_x_ap, element_y_ap;
    mpz_init(scalar_ap);
    mpz_init(element_x_ap);
    mpz_init(element_y_ap);

    std::istringstream stream(std::string(reinterpret_cast<char*>(peer_scalarElement)));
    std::string received_scalar, received_elementX, received_elementY;
    while (std::getline(stream, received_scalar) && std::getline(stream, received_elementX) && std::getline(stream, received_elementY)) {
        mpz_set_str(scalar_ap, received_scalar.c_str(), 10);
        mpz_set_str(element_x_ap, received_elementX.c_str(), 10);
        mpz_set_str(element_y_ap, received_elementY.c_str(), 10);
    }

    // Create a Point object using the parsed coordinates
    Point element_ap(element_x_ap, element_y_ap);

    gmp_printf("---------------------------------------------------\n");
    gmp_printf("Computing shared secret...\n\n");

    cpucycles_reset();
    cpucycles_start();
    std::string sta_token = sta.compute_shared_secret(element_ap, scalar_ap, mac_string_ap);
    cpucycles_stop();
    unsigned int sta_ss_cycles = cpucycles_result();

    mpz_clear(element_x_ap);
    mpz_clear(element_y_ap);
    
    gmp_printf("---------------------------------------------------\n");
    gmp_printf("Confirm Exchange...\n");

    bytes_sent = send(client_fd, sta_token.c_str(), 32, 0);
    clock_gettime(CLOCK_REALTIME, &begin_sending_token);
    std::cout << "\n[STA] Transmitted Data Size " << bytes_sent << " Bytes." << std::endl;
    std::cout << "[STA] STA token sent." << std::endl;

    uint8_t ap_token[32];
    valread = recv(client_fd, ap_token, 32, 0);
    clock_gettime(CLOCK_REALTIME, &end_receiving_token);
    std::cout << "\n[STA] Received Data Size " << valread << " Bytes." << std::endl;
    std::cout << "[STA] AP token received." << std::endl;

    double RTT_TOKEN = (end_receiving_token.tv_sec - begin_sending_token.tv_sec) + (end_receiving_token.tv_nsec - begin_sending_token.tv_nsec) / 1000000000.0 * 1000.0;

    std::string ap_token_string(reinterpret_cast<char*>(ap_token), 32);

    cpucycles_reset();
    cpucycles_start();
    std::string pmk = sta.confirm_exchange(ap_token_string);
    cpucycles_stop();
    unsigned int sta_confirm_cycles = cpucycles_result();

    // write pmk key to file
    ofstream pmkfile("pmk.key");
    pmkfile << pmk;
    pmkfile.close();

    close(client_fd);

    clock_gettime(CLOCK_REALTIME, &end_dragonfly_wall);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_dragonfly_cpu);

    double dragonfly_wall = (end_dragonfly_wall.tv_sec - begin_dragonfly_wall.tv_sec) + (end_dragonfly_wall.tv_nsec - begin_dragonfly_wall.tv_nsec) / 1000000000.0 * 1000.0;
    double dragonfly_cpu = (end_dragonfly_cpu.tv_sec - begin_dragonfly_cpu.tv_sec) + (end_dragonfly_cpu.tv_nsec - begin_dragonfly_cpu.tv_nsec) / 1000000000.0 * 1000.0;

    std::cout << "\nRound Trip Time for MAC Address: " << RTT_MAC << std::endl;
    std::cout << "Round Trip Time for Commit Exchange: " << RTT_SCALAR << std::endl;
    std::cout << "Round Trip Time for Confirm Exchange: " << RTT_TOKEN << std::endl;

    std::cout << "\nTotal CPU Time: " << dragonfly_cpu << " ms" << std::endl;
    std::cout << "Total WALL Time: " << dragonfly_wall << " ms" << std::endl;

    std::cout << "\nDragonfly Initaite: " << sta_initiate_cycles << " CPU Cycles" << std::endl;
    std::cout << "Dragonfly Commit: " << sta_commit_cycles << " CPU Cycles" << std::endl;
    std::cout << "Compute Shared Secret: " << sta_ss_cycles << " CPU Cycles" << std::endl;
    std::cout << "Dragonfly Confirm: " << sta_confirm_cycles << " CPU Cycles" << std::endl;
    std::cout << "Total: " << sta_initiate_cycles + sta_commit_cycles + sta_ss_cycles + sta_confirm_cycles << " CPU Cycles" << std::endl;
    /*
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(c);
    mpz_clear(d);
    */
    return 0;
}
